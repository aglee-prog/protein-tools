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

const SNAPSHOT_KEY: &str = "go_bp_human_proteome_v1";
const BACKGROUND_QUERY: &str = "proteome:UP000005640 AND organism_id:9606";
const FDR_THRESHOLD: f64 = 0.05;
type Error = (StatusCode, String);
fn bad(message: impl Into<String>) -> Error {
    (StatusCode::BAD_REQUEST, message.into())
}
fn upstream(message: impl Into<String>) -> Error {
    (StatusCode::BAD_GATEWAY, message.into())
}
fn internal(message: String) -> Error {
    (StatusCode::INTERNAL_SERVER_ERROR, message)
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
    /// Only human_proteome is supported; omitted means human_proteome.
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
#[derive(Clone, Deserialize, Serialize)]
struct Background {
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
    fn add_page(&mut self, body: &str) -> Result<(), String> {
        let mut lines = body.lines();
        if lines.next() != Some("Entry\tGene Ontology (biological process)") {
            return Err("invalid UniProt GO BP header".into());
        }
        let mut count = 0;
        for line in lines {
            let (accession, annotations) =
                line.split_once('\t').ok_or("invalid UniProt GO BP row")?;
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
    async fn human_go_background(&self) -> Result<Background, String> {
        let endpoint = format!("{}/uniprotkb/search", self.uniprot);
        let base = reqwest::Url::parse(&endpoint).map_err(|_| "invalid UniProt URL")?;
        let mut cursor = None;
        let mut seen = BTreeSet::new();
        let mut expected = None;
        let mut background = Background {
            release: String::new(),
            proteins: BTreeSet::new(),
            terms: BTreeMap::new(),
        };
        for _ in 0..1000 {
            let mut params = vec![
                ("query", BACKGROUND_QUERY.into()),
                ("format", "tsv".into()),
                ("fields", "accession,go_p".into()),
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
    if proteins.is_empty() || !proteins.is_subset(&background.proteins) {
        return Err("input must be a nonempty subset of the human proteome".into());
    }
    let n = proteins.len();
    let population = background.proteins.len();
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
    Ok(EnrichmentResult { method: "overrepresentation".into(), source: "go_bp".into(), background: "human_proteome".into(), background_query: BACKGROUND_QUERY.into(), annotation_policy: "UniProt go_p positive BP associations; all evidence; no additional ancestor propagation; GO root excluded; unannotated proteins included".into(), uniprot_release: background.release.clone(), input_count: n, background_count: population, tested_terms: items.len(), proteins: proteins.into_iter().collect(), items })
}

impl App {
    pub async fn enrich(&self, request: EnrichmentRequest) -> Result<EnrichmentResponse, Error> {
        if request.source.as_deref().is_some_and(|s| s != "go_bp")
            || request
                .background
                .as_deref()
                .is_some_and(|s| s != "human_proteome")
        {
            return Err(bad(
                "only source=go_bp and background=human_proteome are supported",
            ));
        }
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
                    .ok_or((
                        StatusCode::NOT_FOUND,
                        "result_id not found or expired".into(),
                    ))?;
                if value.get("method").and_then(|v| v.as_str()) == Some("overrepresentation") {
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
                    items.into_iter().map(|p| {
                        if !p.found || p.organism.as_deref() != Some("Homo sapiens") {
                            return Err(bad("cached protein result contains unresolved or nonhuman proteins"));
                        }
                        p.uniprot_id.ok_or_else(|| bad("cached protein result lacks an accession"))
                    }).collect::<Result<Vec<_>, _>>()?
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
        let queries: BTreeSet<_> = proteins.iter().map(|p| normalize(p, None).0).collect();
        // Reuse the protein cache and exact resolver without fetching redundant
        // per-input QuickGO annotations: analysis uses one consistent GO snapshot.
        let resolved: Vec<_> = stream::iter(queries.into_iter().map(|query| async move {
            let _permit = self
                .permits
                .acquire()
                .await
                .map_err(|_| internal("server shutting down".into()))?;
            let key = normalize(&query, Some("Homo sapiens"));
            let protein = match self.cache.get(key.clone()).await.map_err(internal)? {
                Some(protein) => Some(protein),
                None => self
                    .upstream
                    .resolve(&key.0, &key.1)
                    .await
                    .map_err(upstream)?,
            }
            .ok_or_else(|| bad(format!("human protein not found: {query}")))?;
            if !protein.found || protein.organism.as_deref() != Some("Homo sapiens") {
                return Err(bad(format!("not a resolved human protein: {query}")));
            }
            protein
                .uniprot_id
                .ok_or_else(|| bad(format!("missing accession: {query}")))
        }))
        .buffered(MAX_CONCURRENT_LOOKUPS)
        .collect()
        .await;
        let proteins: BTreeSet<_> = resolved.into_iter().collect::<Result<_, _>>()?;
        let background = match self
            .cache
            .get_analysis::<Background>(SNAPSHOT_KEY)
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
                let background = self
                    .upstream
                    .human_go_background()
                    .await
                    .map_err(upstream)?;
                self.cache
                    .put_analysis(Some(SNAPSHOT_KEY.into()), background.clone())
                    .await
                    .map_err(internal)?;
                background
            }
        };
        let outside: Vec<_> = proteins.difference(&background.proteins).cloned().collect();
        if !outside.is_empty() {
            return Err(bad(format!(
                "proteins outside human_proteome: {}",
                outside.join(", ")
            )));
        }
        let result = tokio::task::spawn_blocking(move || calculate(&background, proteins))
            .await
            .map_err(|e| internal(e.to_string()))?
            .map_err(internal)?;
        let response = EnrichmentResponse {
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
            release: "test".into(),
            proteins: BTreeSet::new(),
            terms: BTreeMap::new(),
        };
        let mut body = "Entry\tGene Ontology (biological process)\n".to_owned();
        for i in 0..20 {
            let annotation = if i < 7 {
                "process A [GO:0000001]"
            } else if i < 10 {
                "process B [GO:0000002]"
            } else {
                ""
            };
            body.push_str(&format!("P{i:05}\t{annotation}\n"));
        }
        background.add_page(&body).unwrap();
        background
    }
    #[test]
    fn hypergeometric_tail_and_bh_include_zero_hits() {
        let result =
            calculate(&background(), (0..5).map(|i| format!("P{i:05}")).collect()).unwrap();
        // P(X >= 5) = C(7,5) C(13,0) / C(20,5) = 21 / 15504.
        let expected = 21.0 / 15504.0;
        assert!((result.items[0].p_value - expected).abs() < 1e-12);
        assert!((result.items[0].fdr - 2.0 * expected).abs() < 1e-12);
        assert_eq!(result.items[1].p_value, 1.0);
        assert_eq!(result.items[1].hit_count, 0);
        assert_eq!(result.background_count, 20); // includes ten unannotated proteins
        assert_eq!(result.items[0].proteins.len(), 5);
        // A tail with more than one outcome: C(7,3)C(13,2) + C(7,4)C(13,1) + C(7,5).
        let result = calculate(
            &background(),
            ["P00000", "P00001", "P00002", "P00010", "P00011"]
                .map(String::from)
                .into(),
        )
        .unwrap();
        assert!((result.items[0].p_value - 3206.0 / 15504.0).abs() < 1e-12);
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
            bg.add_page("Entry\tGene Ontology (biological process)\nP00000\t\n")
                .is_err()
        );
        assert!(
            bg.add_page("Entry\tGene Ontology (biological process)\nP99999\tbad term\n")
                .is_err()
        );
        assert!(calculate(&bg, BTreeSet::new()).is_err());
        assert!(calculate(&bg, ["ABSENT".into()].into()).is_err());
        let bg = background();
        let all = calculate(&bg, bg.proteins.clone()).unwrap();
        assert!(all.items.iter().all(|t| (t.p_value - 1.0).abs() < 1e-12));
        let none = calculate(&bg, ["P00019".into()].into()).unwrap();
        assert!(none.items.iter().all(|t| t.p_value == 1.0));
        let mut bg = Background {
            release: "test".into(),
            proteins: BTreeSet::new(),
            terms: BTreeMap::new(),
        };
        bg.add_page("Entry\tGene Ontology (biological process)\nP00001\troot [GO:0008150]; test [GO:0000001]; test [GO:0000001]\n").unwrap();
        assert_eq!(bg.terms.len(), 1);
        assert_eq!(bg.terms["GO:0000001"].proteins.len(), 1);
    }
}
