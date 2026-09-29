use crate::{
    api::{App, MAX_CONCURRENT_LOOKUPS, MAX_REQUEST_PROTEINS},
    model::normalize,
    upstream::Upstream,
};
use axum::http::StatusCode;
use futures::{StreamExt, stream};
use serde::{Deserialize, Serialize};
use statrs::distribution::{DiscreteCDF, Hypergeometric};
use std::collections::{BTreeMap, BTreeSet};
use utoipa::ToSchema;

fn background_query(taxon: u64, proteome: &str) -> String {
    format!("proteome:{proteome} AND organism_id:{taxon}")
}
const FDR_THRESHOLD: f64 = 0.05;
pub type Error = (StatusCode, axum::Json<serde_json::Value>);
fn failure(
    status: StatusCode,
    code: &str,
    message: impl Into<String>,
    identifiers: Vec<String>,
) -> Error {
    (
        status,
        axum::Json(
            serde_json::json!({"error": {"code": code, "message": message.into(), "identifiers": identifiers}}),
        ),
    )
}
fn bad(message: impl Into<String>) -> Error {
    failure(StatusCode::BAD_REQUEST, "invalid_request", message, vec![])
}
fn upstream(message: impl Into<String>) -> Error {
    failure(StatusCode::BAD_GATEWAY, "upstream_error", message, vec![])
}
fn internal(message: String) -> Error {
    failure(
        StatusCode::INTERNAL_SERVER_ERROR,
        "internal_error",
        message,
        vec![],
    )
}

#[derive(Deserialize, ToSchema)]
#[serde(deny_unknown_fields)]
pub struct EnrichmentRequest {
    /// Supply exactly one of proteins or result_id. Gene symbols or primary UniProt accessions.
    #[schema(min_items = 1, max_items = 500)]
    pub proteins: Option<Vec<String>>,
    /// ID returned by enrich_proteins or the X-Result-Id header of protein-info.
    pub result_id: Option<String>,
    /// Only go_bp is supported; omitted means go_bp.
    pub source: Option<String>,
    /// NCBI taxonomy ID. Required for a new analysis; inherited from cached enrichment results.
    pub organism_taxon: Option<u64>,
    /// Only proteome is supported; omitted means proteome.
    pub background: Option<String>,
}
#[derive(Clone, Deserialize, Serialize, ToSchema)]
pub struct EnrichmentTerm {
    pub id: String,
    pub name: String,
    pub hit_count: usize,
    pub input_count: usize,
    pub background_hit_count: usize,
    pub background_count: usize,
    pub p_value: f64,
    pub fdr: f64,
    pub proteins: Vec<String>,
}
#[derive(Deserialize, Serialize, ToSchema)]
pub struct EnrichmentResult {
    pub organism_taxon: u64,
    pub input_identifiers: Vec<String>,
    pub identifier_mapping: Vec<IdentifierMapping>,
    pub universe: UniverseMetadata,
    pub method: String,
    pub source: String,
    pub background: String,
    pub background_query: String,
    pub annotation_policy: String,
    pub uniprot_release: String,
    pub input_count: usize,
    pub background_count: usize,
    pub tested_terms: usize,
    pub proteins: Vec<String>,
    pub items: Vec<EnrichmentTerm>,
}
#[derive(Clone, Deserialize, Serialize, ToSchema)]
pub struct IdentifierMapping {
    pub input: String,
    pub canonical: String,
}
#[derive(Clone, Deserialize, Serialize, ToSchema)]
pub struct UniverseMetadata {
    pub proteome_id: String,
    pub unit: String,
    pub selection: String,
    pub proteome_entry_count: usize,
    pub reviewed_entry_count: usize,
    pub unreviewed_entry_count: usize,
    /// Distinct NCBI GeneID cross-references across all proteome entries; not a complete gene census.
    pub unique_gene_id_count: usize,
    pub entries_without_gene_id: usize,
    pub annotated_entry_count: usize,
    pub unannotated_entry_count: usize,
    pub eligible_accessions: Vec<String>,
}
#[derive(Serialize, ToSchema)]
pub struct TopTerm {
    pub id: String,
    pub name: String,
    pub hit_count: usize,
    pub p_value: f64,
    pub fdr: f64,
}
#[derive(Serialize, ToSchema)]
pub struct EnrichmentResponse {
    pub organism_taxon: u64,
    pub result_id: String,
    pub method: String,
    pub source: String,
    pub background: String,
    pub input_count: usize,
    pub background_count: usize,
    pub tested_terms: usize,
    pub significant_terms: usize,
    pub fdr_threshold: f64,
    /// At most ten terms, sorted by FDR, p-value, then GO ID; may be nonsignificant.
    pub top_terms: Vec<TopTerm>,
}
pub const MAX_ENRICHMENT_TERMS: usize = 50;
pub const MAX_DIRECT_ENRICHMENT_BYTES: usize = 32 * 1024;

pub fn is_enrichment(value: &serde_json::Value) -> bool {
    value.get("source").and_then(|v| v.as_str()) == Some("go_bp")
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TermPageQuery {
    #[serde(default = "default_term_limit")]
    pub limit: usize,
    #[serde(default)]
    pub offset: usize,
}
fn default_term_limit() -> usize {
    20
}

#[derive(Serialize, ToSchema)]
pub struct TermStatistics {
    pub id: String,
    pub name: String,
    pub hit_count: usize,
    pub input_count: usize,
    pub background_hit_count: usize,
    pub background_count: usize,
    pub p_value: f64,
    pub fdr: f64,
}
impl From<&EnrichmentTerm> for TermStatistics {
    fn from(t: &EnrichmentTerm) -> Self {
        Self {
            id: t.id.clone(),
            name: t.name.clone(),
            hit_count: t.hit_count,
            input_count: t.input_count,
            background_hit_count: t.background_hit_count,
            background_count: t.background_count,
            p_value: t.p_value,
            fdr: t.fdr,
        }
    }
}
#[derive(Serialize, ToSchema)]
pub struct EnrichmentTermPage {
    pub result_id: String,
    pub total: usize,
    pub limit: usize,
    pub offset: usize,
    pub has_more: bool,
    pub items: Vec<TermStatistics>,
}
#[derive(Serialize, ToSchema)]
pub struct EnrichmentTermDetail {
    pub result_id: String,
    pub organism_taxon: u64,
    pub method: String,
    pub source: String,
    pub background: String,
    pub uniprot_release: String,
    pub term: TermStatistics,
}
#[derive(Serialize, ToSchema)]
pub struct EnrichmentTermProteins {
    pub result_id: String,
    pub term_id: String,
    /// Distinct canonical input hits; aliases do not increase this count.
    pub hit_count: usize,
    pub proteins: Vec<String>,
    /// Original input order, including aliases and repeated inputs.
    pub identifier_mapping: Vec<IdentifierMapping>,
}
impl App {
    pub async fn enrichment_result(&self, id: &str) -> Result<EnrichmentResult, Error> {
        if !valid_result_id(id) {
            return Err(bad("invalid result_id"));
        }
        let value = self
            .cache
            .get_analysis::<serde_json::Value>(id)
            .await
            .map_err(internal)?
            .ok_or_else(|| {
                failure(
                    StatusCode::NOT_FOUND,
                    "result_not_found",
                    "result_id not found or expired",
                    vec![],
                )
            })?;
        if !is_enrichment(&value) {
            return Err(failure(
                StatusCode::BAD_REQUEST,
                "wrong_result_type",
                "result_id is not a GO BP enrichment result",
                vec![],
            ));
        }
        serde_json::from_value(value).map_err(|e| internal(e.to_string()))
    }
}
impl EnrichmentResult {
    pub fn selected_term(&self, id: &str) -> Result<&EnrichmentTerm, Error> {
        if !id
            .strip_prefix("GO:")
            .is_some_and(|s| s.len() == 7 && s.bytes().all(|b| b.is_ascii_digit()))
        {
            return Err(bad(
                "term_id must be a GO ID (GO: followed by seven digits)",
            ));
        }
        self.items.iter().find(|t| t.id == id).ok_or_else(|| {
            failure(
                StatusCode::NOT_FOUND,
                "term_not_found",
                "term_id not found in this result",
                vec![],
            )
        })
    }
    pub fn term_page(
        &mut self,
        id: String,
        query: TermPageQuery,
    ) -> Result<EnrichmentTermPage, Error> {
        if query.limit == 0 {
            return Err(bad("limit must be at least 1"));
        }
        let limit = query.limit.min(MAX_ENRICHMENT_TERMS);
        self.items.sort_by(|a, b| {
            a.fdr
                .total_cmp(&b.fdr)
                .then(a.p_value.total_cmp(&b.p_value))
                .then(a.id.cmp(&b.id))
        });
        let items: Vec<_> = self
            .items
            .iter()
            .skip(query.offset)
            .take(limit)
            .map(TermStatistics::from)
            .collect();
        Ok(EnrichmentTermPage {
            result_id: id,
            total: self.items.len(),
            limit,
            offset: query.offset,
            has_more: query.offset.saturating_add(items.len()) < self.items.len(),
            items,
        })
    }
}

#[derive(Clone, Deserialize, Serialize)]
struct Background {
    proteome_id: String,
    taxon: u64,
    reviewed: usize,
    gene_ids: BTreeSet<String>,
    entries_without_gene_id: usize,
    release: String,
    proteins: BTreeSet<String>,
    terms: BTreeMap<String, BackgroundTerm>,
}
#[derive(Clone, Deserialize, Serialize)]
struct BackgroundTerm {
    name: String,
    proteins: BTreeSet<String>,
}

impl Background {
    fn eligible(&self) -> BTreeSet<String> {
        self.terms
            .values()
            .flat_map(|t| t.proteins.iter().cloned())
            .collect()
    }
    fn add_page(&mut self, body: &str) -> Result<(), String> {
        let mut lines = body.lines();
        if lines.next() != Some("Entry\tReviewed\tGeneID\tGene Ontology (biological process)") {
            return Err("invalid UniProt GO BP header".into());
        }
        let mut count = 0;
        for line in lines {
            let columns: Vec<_> = line.split('\t').collect();
            let [accession, reviewed, genes, annotations] = columns.as_slice() else {
                return Err("invalid UniProt GO BP row".into());
            };
            let (accession, annotations) = (*accession, *annotations);
            match *reviewed {
                "reviewed" => self.reviewed += 1,
                "unreviewed" => (),
                _ => return Err("invalid review status".into()),
            }
            if genes.is_empty() {
                self.entries_without_gene_id += 1;
            }
            self.gene_ids.extend(
                genes
                    .split(';')
                    .map(str::trim)
                    .filter(|g| !g.is_empty())
                    .map(String::from),
            );
            if accession.is_empty()
                || !accession.bytes().all(|c| c.is_ascii_alphanumeric())
                || !self.proteins.insert(accession.to_owned())
            {
                return Err("invalid or duplicate background accession".into());
            }
            count += 1;
            if annotations.is_empty() {
                continue;
            }
            for annotation in annotations.split("; ") {
                let (name, id) = annotation
                    .rsplit_once(" [")
                    .ok_or("invalid GO BP annotation")?;
                let id = id.strip_suffix(']').ok_or("invalid GO BP annotation")?;
                if name.is_empty()
                    || id.len() != 10
                    || !id.starts_with("GO:")
                    || !id[3..].bytes().all(|c| c.is_ascii_digit())
                {
                    return Err("invalid GO BP term".into());
                }
                // The aspect root carries no biological information.
                if id == "GO:0008150" {
                    continue;
                }
                let term = self
                    .terms
                    .entry(id.into())
                    .or_insert_with(|| BackgroundTerm {
                        name: name.into(),
                        proteins: BTreeSet::new(),
                    });
                if term.name != name {
                    return Err("inconsistent GO BP name".into());
                }
                term.proteins.insert(accession.into());
            }
        }
        if count == 0 {
            return Err("empty background page".into());
        }
        Ok(())
    }
}
impl Upstream {
    async fn go_background(&self, taxon: u64) -> Result<Background, String> {
        // Taxonomy searches can include descendants: verify the exact taxon locally.
        let response = self
            .get(
                &format!("{}/proteomes/search", self.uniprot),
                &[
                    ("query", format!("taxonomy_id:{taxon}")),
                    ("format", "json".into()),
                    ("size", "500".into()),
                ],
            )
            .await?;
        if !response.status().is_success() {
            return Err(format!(
                "UniProt proteome discovery returned HTTP {}",
                response.status().as_u16()
            ));
        }
        if response
            .headers()
            .get("link")
            .and_then(|v| v.to_str().ok())
            .is_some_and(|v| v.contains("rel=\"next\""))
        {
            return Err(
                "proteome discovery incomplete; cannot select a unique reference proteome".into(),
            );
        }
        #[derive(Deserialize)]
        struct ProteomeSearch {
            results: Vec<Proteome>,
        }
        #[derive(Deserialize)]
        #[serde(rename_all = "camelCase")]
        struct Proteome {
            id: String,
            proteome_type: String,
            taxonomy: Taxonomy,
        }
        #[derive(Deserialize)]
        #[serde(rename_all = "camelCase")]
        struct Taxonomy {
            taxon_id: u64,
        }
        let search: ProteomeSearch = response
            .json()
            .await
            .map_err(|_| "invalid UniProt proteome discovery response")?;
        let ids: BTreeSet<_> = search
            .results
            .into_iter()
            .filter(|p| p.taxonomy.taxon_id == taxon && p.proteome_type == "Reference proteome")
            .map(|p| p.id)
            .collect();
        if ids.len() != 1 {
            return Err(format!(
                "organism_taxon {taxon} must have exactly one reference proteome; found {}",
                ids.len()
            ));
        }
        let proteome_id = ids.into_iter().next().unwrap();
        if proteome_id.len() != 11
            || !proteome_id.starts_with("UP")
            || !proteome_id[2..].bytes().all(|b| b.is_ascii_digit())
        {
            return Err("invalid reference proteome ID".into());
        }
        let endpoint = format!("{}/uniprotkb/search", self.uniprot);
        let base = reqwest::Url::parse(&endpoint).map_err(|_| "invalid UniProt URL")?;
        let mut cursor = None;
        let mut seen = BTreeSet::new();
        let mut expected = None;
        let mut background = Background {
            proteome_id: proteome_id.clone(),
            taxon,
            reviewed: 0,
            gene_ids: BTreeSet::new(),
            entries_without_gene_id: 0,
            release: String::new(),
            proteins: BTreeSet::new(),
            terms: BTreeMap::new(),
        };
        for _ in 0..1000 {
            let mut params = vec![
                ("query", background_query(taxon, &proteome_id)),
                ("format", "tsv".into()),
                ("fields", "accession,reviewed,xref_geneid,go_p".into()),
                ("size", "500".into()),
            ];
            if let Some(cursor) = cursor.take() {
                params.push(("cursor", cursor));
            }
            let response = self.get(&endpoint, &params).await?;
            if !response.status().is_success() {
                return Err(format!(
                    "UniProt background returned HTTP {}",
                    response.status().as_u16()
                ));
            }
            let headers = response.headers();
            let total = headers
                .get("x-total-results")
                .and_then(|h| h.to_str().ok())
                .and_then(|s| s.parse::<usize>().ok())
                .filter(|n| *n > 0)
                .ok_or("missing background population count")?;
            let release = headers
                .get("x-uniprot-release")
                .and_then(|h| h.to_str().ok())
                .filter(|r| !r.is_empty())
                .ok_or("missing UniProt release")?;
            if expected.is_some_and(|n| n != total)
                || (!background.release.is_empty() && background.release != release)
            {
                return Err("UniProt background changed during pagination".into());
            }
            expected = Some(total);
            background.release = release.into();
            if let Some(link) = headers.get("link") {
                let link = link.to_str().map_err(|_| "invalid background pagination")?;
                // Commas may occur inside the URL (e.g. fields=accession,go_p).
                for part in link.split('<').skip(1) {
                    let (url, attributes) = part
                        .split_once('>')
                        .ok_or("invalid background pagination")?;
                    if !attributes.contains("rel=\"next\"") {
                        continue;
                    }
                    let url = reqwest::Url::parse(url)
                        .map_err(|_| "invalid background pagination URL")?;
                    if url.origin() != base.origin() || url.path() != base.path() {
                        return Err("unexpected background pagination URL".into());
                    }
                    let next = url
                        .query_pairs()
                        .find(|(k, _)| k == "cursor")
                        .map(|(_, v)| v.into_owned())
                        .filter(|s| !s.is_empty())
                        .ok_or("missing background cursor")?;
                    if !seen.insert(next.clone()) {
                        return Err("repeated background cursor".into());
                    }
                    cursor = Some(next);
                }
            }
            background.add_page(
                &response
                    .text()
                    .await
                    .map_err(|_| "invalid UniProt background body")?,
            )?;
            if background.proteins.len() > total {
                return Err("background count exceeds population".into());
            }
            if cursor.is_none() {
                if background.proteins.len() != total || background.terms.is_empty() {
                    return Err("incomplete or unannotated UniProt background".into());
                }
                return Ok(background);
            }
        }
        Err("background pagination limit exceeded".into())
    }
}

fn benjamini_hochberg(items: &mut [EnrichmentTerm]) {
    items.sort_by(|a, b| {
        a.p_value
            .total_cmp(&b.p_value)
            .then_with(|| a.id.cmp(&b.id))
    });
    let count = items.len() as f64;
    let mut adjusted: f64 = 1.0;
    for (index, item) in items.iter_mut().enumerate().rev() {
        adjusted = adjusted.min(item.p_value * count / (index + 1) as f64);
        item.fdr = adjusted;
    }
    items.sort_by(|a, b| {
        a.fdr
            .total_cmp(&b.fdr)
            .then_with(|| a.p_value.total_cmp(&b.p_value))
            .then_with(|| a.id.cmp(&b.id))
    });
}
fn calculate(
    background: &Background,
    proteins: BTreeSet<String>,
) -> Result<EnrichmentResult, String> {
    let eligible = background.eligible();
    if proteins.is_empty() || !proteins.is_subset(&eligible) {
        return Err(
            "input must be a nonempty subset of the GO BP annotated proteome universe".into(),
        );
    }
    let n = proteins.len();
    let population = eligible.len();
    let mut items = Vec::with_capacity(background.terms.len());
    for (id, term) in &background.terms {
        let hits: Vec<_> = proteins.intersection(&term.proteins).cloned().collect();
        let distribution =
            Hypergeometric::new(population as u64, term.proteins.len() as u64, n as u64)
                .map_err(|e| e.to_string())?;
        let p_value = if hits.is_empty() {
            1.0
        } else {
            distribution.sf(hits.len() as u64 - 1).clamp(0.0, 1.0)
        };
        if !p_value.is_finite() {
            return Err("nonfinite enrichment probability".into());
        }
        items.push(EnrichmentTerm {
            id: id.clone(),
            name: term.name.clone(),
            hit_count: hits.len(),
            input_count: n,
            background_hit_count: term.proteins.len(),
            background_count: population,
            p_value,
            fdr: 1.0,
            proteins: hits,
        });
    }
    // Test every annotated background term, including zero-hit terms (p = 1).
    benjamini_hochberg(&mut items);
    Ok(EnrichmentResult {
        organism_taxon: background.taxon,
        input_identifiers: vec![], identifier_mapping: vec![],
        universe: UniverseMetadata {
            proteome_id: background.proteome_id.clone(),
            unit: "UniProtKB primary accession (canonical entry), not gene or isoform".into(),
            selection: "Unique reference proteome for the exact taxon; reviewed and unreviewed; at least one non-root GO BP association; no gene-level collapsing".into(),
            proteome_entry_count: background.proteins.len(),
            reviewed_entry_count: background.reviewed,
            unreviewed_entry_count: background.proteins.len() - background.reviewed,
            unique_gene_id_count: background.gene_ids.len(),
            entries_without_gene_id: background.entries_without_gene_id,
            annotated_entry_count: population,
            unannotated_entry_count: background.proteins.len() - population,
            eligible_accessions: eligible.into_iter().collect(),
        },
        method: "overrepresentation".into(), source: "go_bp".into(), background: "proteome".into(),
        background_query: background_query(background.taxon, &background.proteome_id),
        annotation_policy: "UniProt go_p positive BP associations; all evidence; no additional ancestor propagation; GO root excluded; unannotated entries ineligible; N=annotated_entry_count".into(),
        uniprot_release: background.release.clone(), input_count: n, background_count: population,
        tested_terms: items.len(), proteins: proteins.into_iter().collect(), items,
    })
}

impl App {
    pub async fn enrich(&self, request: EnrichmentRequest) -> Result<EnrichmentResponse, Error> {
        if request.source.as_deref().is_some_and(|s| s != "go_bp")
            || request
                .background
                .as_deref()
                .is_some_and(|s| s != "proteome")
        {
            return Err(bad(
                "only source=go_bp and background=proteome are supported",
            ));
        }
        let mut taxon = request.organism_taxon;
        let mut saved_mapping: Option<Vec<IdentifierMapping>> = None;
        let proteins = match (request.proteins, request.result_id) {
            (Some(proteins), None) => proteins,
            (None, Some(id)) => {
                if !valid_result_id(&id) {
                    return Err(bad("invalid result_id"));
                }
                let value: serde_json::Value = self
                    .cache
                    .get_analysis(&id)
                    .await
                    .map_err(internal)?
                    .ok_or_else(|| {
                        failure(
                            StatusCode::NOT_FOUND,
                            "result_not_found",
                            "result_id not found or expired",
                            vec![],
                        )
                    })?;
                if value.get("method").and_then(|v| v.as_str()) == Some("overrepresentation") {
                    if let Some(saved_taxon) = value.get("organism_taxon").and_then(|v| v.as_u64())
                    {
                        if taxon.is_some_and(|id| id != saved_taxon) {
                            return Err(bad("organism_taxon conflicts with cached result"));
                        }
                        taxon = Some(saved_taxon);
                    }
                    saved_mapping = value
                        .get("identifier_mapping")
                        .cloned()
                        .map(serde_json::from_value)
                        .transpose()
                        .map_err(|_| bad("invalid cached identifier mapping"))?;
                    serde_json::from_value(
                        value
                            .get("proteins")
                            .cloned()
                            .ok_or_else(|| bad("cached result has no protein set"))?,
                    )
                    .map_err(|_| bad("cached result has no protein set"))?
                } else {
                    let items: Vec<crate::model::ProteinResult> = serde_json::from_value(
                        value
                            .get("items")
                            .cloned()
                            .ok_or_else(|| bad("cached result has no protein records"))?,
                    )
                    .map_err(|_| bad("cached result has no protein records"))?;
                    let invalid: Vec<_> = items
                        .iter()
                        .filter(|p| !p.found || p.uniprot_id.is_none())
                        .map(|p| p.query.clone())
                        .collect();
                    if !invalid.is_empty() {
                        return Err(failure(
                            StatusCode::BAD_REQUEST,
                            "resolution_failed",
                            "cached result contains unresolved proteins",
                            invalid,
                        ));
                    }
                    let mapping: Vec<_> = items
                        .into_iter()
                        .map(|p| IdentifierMapping {
                            input: p.query,
                            canonical: p.uniprot_id.unwrap(),
                        })
                        .collect();
                    let proteins = mapping.iter().map(|p| p.canonical.clone()).collect();
                    saved_mapping = Some(mapping);
                    proteins
                }
            }
            _ => return Err(bad("supply exactly one of proteins or result_id")),
        };
        if proteins.is_empty()
            || proteins.len() > MAX_REQUEST_PROTEINS
            || proteins
                .iter()
                .any(|p| p.trim().is_empty() || p.len() > 128 || p.chars().any(char::is_control))
        {
            return Err(bad("proteins must contain 1 to 500 valid identifiers"));
        }
        let taxon = taxon.filter(|id| *id > 0).ok_or_else(|| failure(
            StatusCode::BAD_REQUEST, "organism_required", "provide organism_taxon (NCBI taxonomy ID) to select the resolution context and proteome universe", proteins.clone()))?;
        let queries: BTreeSet<_> = proteins.iter().map(|p| normalize(p, None).0).collect();
        // Taxon-scoped cache namespace avoids stale name/alias resolution from v1.
        // Both enrichment and lookup use the same resolver and ranking policy.
        let resolved: Vec<_> = stream::iter(queries.into_iter().map(|query| async move {
            let outcome = async {
                let _permit = self
                    .permits
                    .acquire()
                    .await
                    .map_err(|_| "server shutting down".to_owned())?;
                let key = (query.clone(), format!("taxon:{taxon}:resolver_v2"));
                let protein = match self.cache.get(key.clone()).await? {
                    Some(protein) => Some(protein),
                    None => {
                        let protein = self.upstream.resolve_taxon(&query, taxon).await?;
                        if let Some(protein) = &protein {
                            self.cache.put(key, protein.clone()).await?;
                        }
                        protein
                    }
                };
                Ok::<_, String>(protein.filter(|p| p.found).and_then(|p| p.uniprot_id))
            }
            .await;
            (query, outcome)
        }))
        .buffered(MAX_CONCURRENT_LOOKUPS)
        .collect()
        .await;
        let mut canonical = BTreeMap::new();
        let mut issues = vec![];
        let mut upstream_failed = false;
        for (identifier, outcome) in resolved {
            match outcome {
                Ok(Some(id)) => {
                    canonical.insert(identifier, id);
                }
                Ok(None) => {
                    issues.push(serde_json::json!({"identifier":identifier,"reason":"unresolved"}))
                }
                Err(message) => {
                    let reason = if message.contains("ambiguous")
                        || message.contains("candidate set too large")
                    {
                        "ambiguous"
                    } else {
                        upstream_failed = true;
                        "upstream_error"
                    };
                    issues.push(serde_json::json!({"identifier":identifier,"reason":reason,"message":message}));
                }
            }
        }
        if !issues.is_empty() {
            return Err((
                if upstream_failed {
                    StatusCode::BAD_GATEWAY
                } else {
                    StatusCode::BAD_REQUEST
                },
                axum::Json(serde_json::json!({"error":{
                    "code":"resolution_failed", "organism_taxon":taxon, "message":"Every identifier must resolve; no analysis performed", "identifiers":issues
                }})),
            ));
        }
        let mapping = saved_mapping.unwrap_or_else(|| {
            proteins
                .iter()
                .map(|p| IdentifierMapping {
                    input: p.clone(),
                    canonical: canonical[&normalize(p, None).0].clone(),
                })
                .collect()
        });
        let proteins: BTreeSet<_> = canonical.into_values().collect();
        let snapshot_key = format!("go_bp_proteome_taxon_{taxon}_v2");
        let background = match self
            .cache
            .get_analysis::<Background>(&snapshot_key)
            .await
            .map_err(internal)?
        {
            Some(background) => background,
            None => {
                let _permit = self
                    .permits
                    .acquire()
                    .await
                    .map_err(|_| internal("server shutting down".into()))?;
                let background = self.upstream.go_background(taxon).await.map_err(upstream)?;
                self.cache
                    .put_analysis(Some(snapshot_key), background.clone())
                    .await
                    .map_err(internal)?;
                background
            }
        };
        let eligible = background.eligible();
        let outside: Vec<_> = proteins.difference(&eligible).cloned().collect();
        if !outside.is_empty() {
            return Err(failure(
                StatusCode::BAD_REQUEST,
                "outside_annotation_universe",
                "proteins outside the selected taxon's proteome or without non-root GO BP annotations; no analysis performed",
                outside,
            ));
        }
        let mut result = tokio::task::spawn_blocking(move || calculate(&background, proteins))
            .await
            .map_err(|e| internal(e.to_string()))?
            .map_err(internal)?;
        result.input_identifiers = mapping.iter().map(|m| m.input.clone()).collect();
        result.identifier_mapping = mapping;
        let response = EnrichmentResponse {
            organism_taxon: taxon,
            result_id: String::new(),
            method: result.method.clone(),
            source: result.source.clone(),
            background: result.background.clone(),
            input_count: result.input_count,
            background_count: result.background_count,
            tested_terms: result.tested_terms,
            significant_terms: result
                .items
                .iter()
                .filter(|t| t.fdr <= FDR_THRESHOLD)
                .count(),
            fdr_threshold: FDR_THRESHOLD,
            top_terms: result
                .items
                .iter()
                .take(10)
                .map(|t| TopTerm {
                    id: t.id.clone(),
                    name: t.name.clone(),
                    hit_count: t.hit_count,
                    p_value: t.p_value,
                    fdr: t.fdr,
                })
                .collect(),
        };
        let result_id = self
            .cache
            .put_analysis(None, result)
            .await
            .map_err(internal)?;
        Ok(EnrichmentResponse {
            result_id,
            ..response
        })
    }
}
pub fn valid_result_id(id: &str) -> bool {
    id.strip_prefix("result_")
        .is_some_and(|s| s.len() == 32 && s.bytes().all(|c| c.is_ascii_hexdigit()))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn background() -> Background {
        let mut background = Background {
            proteome_id: "UP000005640".into(),
            taxon: 9606,
            reviewed: 0,
            gene_ids: BTreeSet::new(),
            entries_without_gene_id: 0,
            release: "test".into(),
            proteins: BTreeSet::new(),
            terms: BTreeMap::new(),
        };
        let mut body = "Entry\tReviewed\tGeneID\tGene Ontology (biological process)\n".to_owned();
        for i in 0..20 {
            let annotation = if i < 7 {
                "process A [GO:0000001]"
            } else if i < 10 {
                "process B [GO:0000002]"
            } else {
                ""
            };
            body.push_str(&format!("P{i:05}\treviewed\t{i};\t{annotation}\n"));
        }
        background.add_page(&body).unwrap();
        background
    }
    #[test]
    fn hypergeometric_tail_and_bh_include_zero_hits() {
        let result =
            calculate(&background(), (0..5).map(|i| format!("P{i:05}")).collect()).unwrap();
        // N=10 annotated entries: P(X >= 5) = C(7,5) / C(10,5).
        let expected = 21.0 / 252.0;
        assert!((result.items[0].p_value - expected).abs() < 1e-12);
        assert!((result.items[0].fdr - 2.0 * expected).abs() < 1e-12);
        assert_eq!(result.items[1].p_value, 1.0);
        assert_eq!(result.items[1].hit_count, 0);
        assert_eq!(result.background_count, 10); // excludes ten unannotated proteins
        assert_eq!(result.items[0].proteins.len(), 5);
        // A tail with more than one outcome: C(7,3)C(3,2) + C(7,4)C(3,1) + C(7,5).
        let result = calculate(
            &background(),
            ["P00000", "P00001", "P00002", "P00007", "P00008"]
                .map(String::from)
                .into(),
        )
        .unwrap();
        assert!(
            (result
                .items
                .iter()
                .find(|t| t.id == "GO:0000001")
                .unwrap()
                .p_value
                - 231.0 / 252.0)
                .abs()
                < 1e-12
        );
    }
    #[test]
    fn bh_monotonicity_ties_and_bounds() {
        let mut terms = calculate(&background(), ["P00000".into()].into())
            .unwrap()
            .items;
        let template = terms[0].clone();
        terms = [0.01, 0.04, 0.03, 0.002, 1.0]
            .into_iter()
            .enumerate()
            .map(|(i, p)| EnrichmentTerm {
                id: i.to_string(),
                p_value: p,
                ..template.clone()
            })
            .collect();
        benjamini_hochberg(&mut terms);
        for (t, expected) in terms.iter().zip([0.01, 0.025, 0.05, 0.05, 1.0]) {
            assert!((t.fdr - expected).abs() < 1e-12);
        }
    }
    #[test]
    fn parsing_and_population_edges() {
        let mut bg = background();
        assert!(
            bg.add_page("Entry\tReviewed\tGeneID\tGene Ontology (biological process)\nP00000\treviewed\t\t\n")
                .is_err()
        );
        assert!(
            bg.add_page("Entry\tReviewed\tGeneID\tGene Ontology (biological process)\nP99999\treviewed\t\tbad term\n")
                .is_err()
        );
        assert!(calculate(&bg, BTreeSet::new()).is_err());
        assert!(calculate(&bg, ["ABSENT".into()].into()).is_err());
        let bg = background();
        let all = calculate(&bg, bg.eligible()).unwrap();
        assert!(all.items.iter().all(|t| (t.p_value - 1.0).abs() < 1e-12));
        assert!(calculate(&bg, ["P00019".into()].into()).is_err());
        let mut bg = Background {
            proteome_id: "UP000005640".into(),
            taxon: 9606,
            reviewed: 0,
            gene_ids: BTreeSet::new(),
            entries_without_gene_id: 0,
            release: "test".into(),
            proteins: BTreeSet::new(),
            terms: BTreeMap::new(),
        };
        bg.add_page("Entry\tReviewed\tGeneID\tGene Ontology (biological process)\nP00001\treviewed\t\troot [GO:0008150]; test [GO:0000001]; test [GO:0000001]\n").unwrap();
        assert_eq!(bg.terms.len(), 1);
        assert_eq!(bg.terms["GO:0000001"].proteins.len(), 1);
    }
}
