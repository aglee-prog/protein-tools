//! Multi-organism KEGG lookup, independent of the HTTP/tool transport.
use crate::{cache::Cache, uniprot::ProteinTaxonomy, upstream::Upstream};
use reqwest::Client;
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    time::Duration,
};
use tokio::{sync::Mutex, time::Instant};
use utoipa::ToSchema;

const BATCH_SIZE: usize = 100;
// Shared by every provider and request in this process, including retries.
static NEXT_REQUEST: Mutex<Option<Instant>> = Mutex::const_new(None);

#[derive(Deserialize, ToSchema)]
pub struct KeggRequest {
    /// Canonical UniProt accessions, 1 to 500, from any mix of organisms. Send the entire set in one call.
    #[schema(min_items = 1, max_items = 500)]
    pub proteins: Vec<String>,
    /// Optional organism-specific pathway ID (the path: prefix is also accepted).
    /// Filters the pathway groups only; protein results retain all available associations.
    pub pathway_id: Option<String>,
    /// Minimum distinct supplied accessions in a pathway group (default 1).
    /// Use 2 for shared pathways, or the distinct input count for intersection.
    #[schema(minimum = 1, maximum = 500)]
    pub min_proteins: Option<usize>,
}
#[derive(Clone, Debug, Serialize, ToSchema)]
pub struct KeggAssociation {
    pub kegg_gene_id: String,
    pub pathway_id: String,
    /// Null when name retrieval failed; the authoritative association is retained.
    pub pathway_name: Option<String>,
}
#[derive(Clone, Debug, Serialize, ToSchema)]
#[serde(rename_all = "snake_case")]
pub enum KeggStatus {
    Mapped,
    NoMapping,
    UnsupportedOrganism,
    NoPathways,
    Error,
}
#[derive(Clone, Debug, Serialize, ToSchema)]
pub struct KeggProtein {
    pub query: String,
    pub uniprot_id: String,
    /// NCBI taxonomy ID reported by UniProt.
    pub ncbi_taxonomy_id: Option<u64>,
    /// Scientific organism name reported by UniProt, when available.
    pub organism_name: Option<String>,
    /// All KEGG organism codes explicitly linked to this exact taxonomy ID.
    pub kegg_organism_codes: Vec<String>,
    pub kegg_gene_ids: Vec<String>,
    pub pathways: Vec<KeggAssociation>,
    pub status: KeggStatus,
    /// True only when all required stages were served from persistent cache.
    pub cached: bool,
    /// Non-null means results for this protein may be incomplete.
    pub error: Option<String>,
    pub source: String,
}
#[derive(Debug, Serialize, ToSchema)]
pub struct KeggPathwayGroup {
    pub pathway_id: String,
    pub pathway_name: Option<String>,
    pub uniprot_ids: Vec<String>,
    pub protein_count: usize,
    pub source: String,
}
#[derive(Debug, Serialize, ToSchema)]
pub struct KeggResponse {
    /// False if any input failed. Groups then reflect only known associations.
    pub complete: bool,
    /// One result per input, in input order, including duplicates and failures.
    pub proteins: Vec<KeggProtein>,
    /// Sorted groups; duplicate inputs and multiple genes do not inflate counts.
    pub pathways: Vec<KeggPathwayGroup>,
}

pub fn accession(query: &str) -> Result<String, String> {
    let value = query.trim().to_ascii_uppercase();
    let b = value.as_bytes();
    let valid = matches!(b.len(), 6 | 10)
        && b.iter().all(u8::is_ascii_alphanumeric)
        && b[0].is_ascii_uppercase()
        && b[1].is_ascii_digit()
        && if b"OPQ".contains(&b[0]) {
            b.len() == 6 && b[5].is_ascii_digit()
        } else {
            b[2..]
                .chunks(4)
                .all(|c| c[0].is_ascii_uppercase() && c[3].is_ascii_digit())
        };
    if valid {
        Ok(value)
    } else {
        Err("expected a canonical UniProt accession (gene symbols and isoform suffixes are not supported)".into())
    }
}
pub fn pathway_id(value: &str) -> Result<String, String> {
    let value = value.trim().to_ascii_lowercase();
    let value = value.strip_prefix("path:").unwrap_or(&value);
    if value.is_ascii() && value.len() > 5 {
        let (organism, number) = value.split_at(value.len() - 5);
        if organism_code(organism)
            && !["map", "ko", "ec", "rn"].contains(&organism)
            && number.bytes().all(|b| b.is_ascii_digit())
        {
            return Ok(value.into());
        }
    }
    Err("expected an organism-specific KEGG pathway ID".into())
}
fn organism_code(value: &str) -> bool {
    (3..=4).contains(&value.len()) && value.bytes().all(|b| b.is_ascii_lowercase())
}
fn gene_id(value: &str, organism: &str) -> bool {
    value.split_once(':').is_some_and(|(prefix, gene)| {
        prefix == organism
            && !gene.is_empty()
            && gene
                .bytes()
                .all(|b| b.is_ascii_alphanumeric() || b"_.-".contains(&b))
    })
}
fn organism_pathway(value: &str, organism: &str) -> Result<String, String> {
    let id = pathway_id(value)?;
    if id.strip_prefix(organism).is_some_and(|n| n.len() == 5) {
        Ok(id)
    } else {
        Err("KEGG returned a pathway for another organism".into())
    }
}

#[derive(Clone, Copy)]
enum Stage {
    Taxonomy,
    Conversion,
    Links,
    Names,
}
impl Stage {
    fn key(self) -> &'static str {
        match self {
            Self::Taxonomy => "taxonomy-codes-v1",
            Self::Conversion => "conversion-v1",
            Self::Links => "links-v1",
            Self::Names => "names-v1",
        }
    }
    fn path(self, organism: &str, ids: &[String]) -> String {
        match self {
            Self::Taxonomy => format!(
                "link/genome/{}",
                ids.iter()
                    .map(|id| format!("taxid:{id}"))
                    .collect::<Vec<_>>()
                    .join("+")
            ),
            Self::Conversion => format!(
                "conv/{organism}/{}",
                ids.iter()
                    .map(|id| format!("up:{id}"))
                    .collect::<Vec<_>>()
                    .join("+")
            ),
            Self::Links => format!("link/pathway/{}", ids.join("+")),
            Self::Names => format!("list/pathway/{organism}"),
        }
    }
}
// Parse an entire response before writing any of it. Unexpected identifiers or
// malformed bodies must never turn into cached negative biological results.
fn parse(
    stage: Stage,
    organism: &str,
    body: &str,
    ids: &[String],
) -> Result<BTreeMap<String, Vec<String>>, String> {
    let mut rows: BTreeMap<String, BTreeSet<String>> =
        ids.iter().map(|id| (id.clone(), BTreeSet::new())).collect();
    for line in body.lines().filter(|s| !s.trim().is_empty()) {
        let (left, right) = line
            .split_once('\t')
            .ok_or("invalid KEGG tabular response")?;
        if right.is_empty() || right.contains('\t') {
            return Err("invalid KEGG tabular response".into());
        }
        let (key, value) = match stage {
            Stage::Taxonomy => {
                let taxid = left
                    .strip_prefix("taxid:")
                    .ok_or("invalid KEGG taxonomy ID")?;
                let code = right
                    .strip_prefix("gn:")
                    .ok_or("invalid KEGG genome identifier")?;
                if !organism_code(code) {
                    return Err("invalid KEGG organism code".into());
                }
                (taxid.to_owned(), code.to_owned())
            }
            Stage::Conversion => {
                let key = left
                    .strip_prefix("up:")
                    .or_else(|| left.strip_prefix("uniprot:"))
                    .ok_or("invalid KEGG UniProt identifier")?;
                if !gene_id(right, organism) {
                    return Err("invalid KEGG gene identifier".into());
                }
                (key.to_owned(), right.to_owned())
            }
            Stage::Links => {
                if !gene_id(left, organism) {
                    return Err("invalid KEGG gene identifier".into());
                }
                (left.to_owned(), organism_pathway(right, organism)?)
            }
            Stage::Names => (organism_pathway(left, organism)?, right.to_owned()),
        };
        if let Some(values) = rows.get_mut(&key) {
            values.insert(value);
        } else if !matches!(stage, Stage::Names) {
            return Err("KEGG returned an unrequested identifier".into());
        }
        // The organism name catalog also contains pathways not in this request.
        // Validate those rows above, but only retain requested memberships.
    }
    Ok(rows
        .into_iter()
        .map(|(id, values)| (id, values.into_iter().collect()))
        .collect())
}

fn valid_name(values: &[String]) -> bool {
    values.len() == 1 && !values[0].trim().is_empty()
}

fn retry_after(value: &str) -> Option<u64> {
    value.parse().ok().or_else(|| {
        httpdate::parse_http_date(value).ok().map(|date| {
            date.duration_since(std::time::SystemTime::now())
                .unwrap_or_default()
                .as_secs()
                .saturating_add(1)
        })
    })
}

#[derive(Clone)]
pub struct Kegg {
    client: Client,
    base_url: String,
}
#[derive(Clone)]
struct Values {
    values: Vec<String>,
    cached: bool,
}
type StageResults = BTreeMap<String, Result<Values, String>>;
impl Kegg {
    pub fn new(base_url: String) -> Result<Self, reqwest::Error> {
        Ok(Self {
            client: Client::builder()
                .connect_timeout(Duration::from_secs(5))
                .timeout(Duration::from_secs(30))
                .redirect(reqwest::redirect::Policy::none())
                .user_agent("protein-tools/0.1.0")
                .build()?,
            base_url,
        })
    }
    async fn get(&self, path: &str) -> Result<String, String> {
        let mut error = "KEGG request failed".to_owned();
        for attempt in 0..3 {
            // Hold through the HTTP response so Retry-After can delay all callers.
            let mut next = NEXT_REQUEST.lock().await;
            if let Some(deadline) = *next {
                if deadline.saturating_duration_since(Instant::now()) > Duration::from_secs(30) {
                    return Err("KEGG cooldown active; retry later".into());
                }
                tokio::time::sleep_until(deadline).await;
            }
            *next = Some(Instant::now() + Duration::from_millis(500));
            let response = self
                .client
                .get(format!("{}/{}", self.base_url.trim_end_matches('/'), path))
                .header("Accept", "text/plain")
                .send()
                .await;
            match response {
                Ok(response) if response.status().is_success() => {
                    return response
                        .text()
                        .await
                        .map_err(|_| "invalid KEGG response body".into());
                }
                Ok(response) => {
                    let status = response.status();
                    error = format!("KEGG returned HTTP {}", status.as_u16());
                    if status.as_u16() != 429 && !status.is_server_error() {
                        return Err(error);
                    }
                    let delay = response
                        .headers()
                        .get("retry-after")
                        .and_then(|v| v.to_str().ok())
                        .and_then(retry_after)
                        .unwrap_or(1 << attempt);
                    // A long server cooldown is preserved without keeping this request asleep indefinitely.
                    *next = Instant::now().checked_add(Duration::from_secs(delay.max(1)));
                    if delay > 30 {
                        return Err(format!("{error}; retry later"));
                    }
                }
                Err(_) => {
                    error = "KEGG network request failed".into();
                    *next = Some(Instant::now() + Duration::from_secs(1 << attempt));
                }
            }
        }
        Err(error)
    }
    async fn stage(
        &self,
        cache: &Cache,
        organism: &str,
        stage: Stage,
        ids: BTreeSet<String>,
    ) -> StageResults {
        let mut results = BTreeMap::new();
        let mut misses = Vec::new();
        for id in ids {
            match cache.get_kegg(organism, stage.key(), &id).await {
                Ok(Some(values)) if !matches!(stage, Stage::Names) || valid_name(&values) => {
                    results.insert(
                        id,
                        Ok(Values {
                            values,
                            cached: true,
                        }),
                    );
                }
                other => {
                    if let Err(error) = other {
                        tracing::warn!(event="kegg_cache_read_failure", %error);
                    }
                    misses.push(id);
                }
            }
        }
        // Fetch the organism catalog once for all missing names. KEGG's
        // explicit-entry list operation does not return organism pathway names.
        let size = if matches!(stage, Stage::Names) {
            misses.len().max(1)
        } else {
            BATCH_SIZE
        };
        for chunk in misses.chunks(size) {
            let parsed = match self.get(&stage.path(organism, chunk)).await {
                Ok(body) => parse(stage, organism, &body, chunk),
                Err(error) => Err(error),
            };
            match parsed {
                Ok(rows) => {
                    for (id, values) in rows {
                        if matches!(stage, Stage::Names) && !valid_name(&values) {
                            results
                                .insert(id, Err("KEGG pathway name missing or ambiguous".into()));
                            continue;
                        }
                        if let Err(error) = cache
                            .put_kegg(organism, stage.key(), &id, values.clone())
                            .await
                        {
                            tracing::warn!(event="kegg_cache_write_failure", %error);
                        }
                        results.insert(
                            id,
                            Ok(Values {
                                values,
                                cached: false,
                            }),
                        );
                    }
                }
                Err(error) => {
                    for id in chunk {
                        results.insert(id.clone(), Err(error.clone()));
                    }
                }
            }
        }
        results
    }
    async fn resolve_taxonomies(
        &self,
        cache: &Cache,
        upstream: &Upstream,
        ids: BTreeSet<String>,
    ) -> BTreeMap<String, Result<(Option<ProteinTaxonomy>, bool), String>> {
        let mut results = BTreeMap::new();
        let mut misses = Vec::new();
        for id in ids {
            match cache
                .get_kegg_value::<Option<ProteinTaxonomy>>("", "uniprot-taxonomy-v1", &id)
                .await
            {
                Ok(Some(value)) => {
                    results.insert(id, Ok((value, true)));
                }
                other => {
                    if let Err(error) = other {
                        tracing::warn!(event="kegg_cache_read_failure", %error);
                    }
                    misses.push(id);
                }
            }
        }
        for chunk in misses.chunks(BATCH_SIZE) {
            match upstream.taxonomies(chunk).await {
                Ok(values) => {
                    for (id, value) in values {
                        if let Err(error) = cache
                            .put_kegg_value("", "uniprot-taxonomy-v1", &id, value.clone())
                            .await
                        {
                            tracing::warn!(event="kegg_cache_write_failure", %error);
                        }
                        results.insert(id, Ok((value, false)));
                    }
                }
                Err(error) => {
                    for id in chunk {
                        results.insert(id.clone(), Err(error.clone()));
                    }
                }
            }
        }
        results
    }
    pub async fn lookup(
        &self,
        cache: &Cache,
        upstream: &Upstream,
        request: KeggRequest,
    ) -> KeggResponse {
        let ids: BTreeSet<String> = request
            .proteins
            .iter()
            .filter_map(|s| accession(s).ok())
            .collect();
        let taxonomies = self.resolve_taxonomies(cache, upstream, ids.clone()).await;
        let taxids = taxonomies
            .values()
            .filter_map(|r| r.as_ref().ok())
            .filter_map(|(t, _)| t.as_ref()?.ncbi_taxonomy_id.map(|id| id.to_string()))
            .collect();
        let codes = self.stage(cache, "", Stage::Taxonomy, taxids).await;
        let mut results = BTreeMap::new();
        let mut groups: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
        for id in ids {
            let mut protein = empty_protein(&id);
            match &taxonomies[&id] {
                Err(error) => fail(&mut protein, error),
                Ok((taxonomy, cached)) => {
                    protein.cached = *cached;
                    if let Some(taxonomy) = taxonomy {
                        protein.ncbi_taxonomy_id = taxonomy.ncbi_taxonomy_id;
                        protein.organism_name = taxonomy.organism_name.clone();
                        protein.status = KeggStatus::UnsupportedOrganism;
                        if let Some(taxid) = taxonomy.ncbi_taxonomy_id {
                            match &codes[&taxid.to_string()] {
                                Err(error) => fail(&mut protein, error),
                                Ok(value) => {
                                    protein.cached &= value.cached;
                                    protein.kegg_organism_codes = value.values.clone();
                                    for code in &value.values {
                                        groups.entry(code.clone()).or_default().insert(id.clone());
                                    }
                                }
                            }
                        }
                    }
                }
            }
            results.insert(id, protein);
        }
        for (organism, ids) in groups {
            let mappings = self.stage(cache, &organism, Stage::Conversion, ids).await;
            let links = self
                .stage(cache, &organism, Stage::Links, successful_values(&mappings))
                .await;
            let names = self
                .stage(cache, &organism, Stage::Names, successful_values(&links))
                .await;
            for (id, mapping) in mappings {
                let protein = results.get_mut(&id).expect("resolved accession");
                match mapping {
                    Err(error) => fail(protein, &error),
                    Ok(mapping) => {
                        protein.cached &= mapping.cached;
                        protein.kegg_gene_ids.extend(mapping.values.clone());
                        for gene in &mapping.values {
                            match &links[gene] {
                                Err(error) => fail(protein, error),
                                Ok(paths) => {
                                    protein.cached &= paths.cached;
                                    for path in &paths.values {
                                        let name = match &names[path] {
                                            Err(error) => {
                                                fail(protein, error);
                                                None
                                            }
                                            Ok(value) => {
                                                protein.cached &= value.cached;
                                                value.values.first().cloned()
                                            }
                                        };
                                        protein.pathways.push(KeggAssociation {
                                            kegg_gene_id: gene.clone(),
                                            pathway_id: path.clone(),
                                            pathway_name: name,
                                        });
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        let proteins: Vec<_> = request
            .proteins
            .into_iter()
            .map(|query| {
                let mut protein = match accession(&query) {
                    Ok(id) => results[&id].clone(),
                    Err(error) => {
                        let mut p = empty_protein(&query);
                        fail(&mut p, &error);
                        p
                    }
                };
                protein.query = query;
                protein.status = if protein.error.is_some() {
                    KeggStatus::Error
                } else if matches!(protein.status, KeggStatus::UnsupportedOrganism)
                    && protein.kegg_organism_codes.is_empty()
                {
                    KeggStatus::UnsupportedOrganism
                } else if protein.kegg_gene_ids.is_empty() {
                    KeggStatus::NoMapping
                } else if protein.pathways.is_empty() {
                    KeggStatus::NoPathways
                } else {
                    KeggStatus::Mapped
                };
                protein
            })
            .collect();
        let pathways = aggregate(
            &proteins,
            request.pathway_id.as_deref(),
            request.min_proteins.unwrap_or(1),
        );
        KeggResponse {
            complete: proteins.iter().all(|p| p.error.is_none()),
            proteins,
            pathways,
        }
    }
}
fn empty_protein(query: &str) -> KeggProtein {
    KeggProtein {
        query: query.into(),
        uniprot_id: query.trim().to_ascii_uppercase(),
        ncbi_taxonomy_id: None,
        organism_name: None,
        kegg_organism_codes: vec![],
        kegg_gene_ids: vec![],
        pathways: vec![],
        status: KeggStatus::NoMapping,
        cached: true,
        error: None,
        source: "KEGG".into(),
    }
}
fn fail(protein: &mut KeggProtein, error: &str) {
    protein.error = Some(error.into());
    protein.cached = false;
}

fn successful_values(results: &StageResults) -> BTreeSet<String> {
    results
        .values()
        .filter_map(|r| r.as_ref().ok())
        .flat_map(|r| r.values.iter().cloned())
        .collect()
}
fn aggregate(
    proteins: &[KeggProtein],
    filter: Option<&str>,
    minimum: usize,
) -> Vec<KeggPathwayGroup> {
    let mut groups: BTreeMap<String, (Option<String>, BTreeSet<String>)> = BTreeMap::new();
    for protein in proteins {
        for path in &protein.pathways {
            if filter.is_some_and(|id| id != path.pathway_id) {
                continue;
            }
            let group = groups
                .entry(path.pathway_id.clone())
                .or_insert_with(|| (path.pathway_name.clone(), BTreeSet::new()));
            group.1.insert(protein.uniprot_id.clone());
        }
    }
    groups
        .into_iter()
        .filter(|(_, (_, ids))| ids.len() >= minimum)
        .map(|(pathway_id, (pathway_name, ids))| KeggPathwayGroup {
            pathway_id,
            pathway_name,
            protein_count: ids.len(),
            uniprot_ids: ids.into_iter().collect(),
            source: "KEGG".into(),
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn identifiers() {
        for id in ["P04637", "Q9Y243", "A0A024RBG1", "B2R8Q1", "O14920"] {
            assert_eq!(accession(&format!(" {} ", id.to_lowercase())).unwrap(), id);
        }
        for id in [
            "",
            "TP53",
            "P04637-2",
            "up:P04637",
            "P04637+P00533",
            "P0463X",
            "A01AAA",
            "é12345",
        ] {
            assert!(accession(id).is_err(), "{id}");
        }
        assert_eq!(pathway_id(" path:HSA04115 ").unwrap(), "hsa04115");
        for id in ["map04115", "mmu:04115", "hsa:04115", "hsa04115+", "hsaé123"] {
            assert!(pathway_id(id).is_err());
        }
    }
    #[test]
    fn tabular_parsing_is_strict_and_preserves_multiple_mappings() {
        let ids = vec!["P04637".into(), "P00533".into()];
        let rows = parse(
            Stage::Conversion,
            "hsa",
            "up:P04637\thsa:7157\r\nup:P04637\thsa:7157\n uniprot:P04637\thsa:1",
            &ids,
        );
        assert!(rows.is_err());
        let rows = parse(
            Stage::Conversion,
            "hsa",
            "up:P04637\thsa:7157\r\nuniprot:P04637\thsa:1\nup:P04637\thsa:7157\n",
            &ids,
        )
        .unwrap();
        assert_eq!(rows["P04637"], ["hsa:1", "hsa:7157"]);
        assert!(rows["P00533"].is_empty());
        for body in [
            "<html>error</html>",
            "up:P04637\tmmu:1",
            "up:Q99999\thsa:1",
            "up:P04637\thsa:1\textra",
        ] {
            assert!(parse(Stage::Conversion, "hsa", body, &ids).is_err());
        }
        let genes = vec!["hsa:7157".into()];
        assert_eq!(
            parse(Stage::Links, "hsa", "hsa:7157\tpath:hsa04115", &genes).unwrap()["hsa:7157"],
            ["hsa04115"]
        );
        assert!(parse(Stage::Links, "hsa", "hsa:7157\tpath:mmu04115", &genes).is_err());
        let paths = vec!["hsa04115".into()];
        assert_eq!(
            parse(
                Stage::Names,
                "hsa",
                "hsa04115\tp53 signaling pathway - Homo sapiens (human)",
                &paths
            )
            .unwrap()["hsa04115"],
            ["p53 signaling pathway - Homo sapiens (human)"]
        );
        assert!(parse(Stage::Names, "hsa", "path:hsa04115\tp53", &paths).is_ok());
    }
    #[test]
    fn retry_after_formats() {
        assert_eq!(retry_after("2"), Some(2));
        assert_eq!(retry_after("invalid"), None);
        let date = httpdate::fmt_http_date(std::time::SystemTime::now() + Duration::from_secs(10));
        assert!((9..=11).contains(&retry_after(&date).unwrap()));
    }
}
