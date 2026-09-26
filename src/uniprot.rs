use crate::{model::ProteinResult, upstream::Upstream};
use serde::Deserialize;

#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Record {
    primary_accession: String,
    entry_type: String,
    protein_description: Description,
    #[serde(default)]
    genes: Vec<Gene>,
    organism: Organism,
    #[serde(default)]
    comments: Vec<Comment>,
}
#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct Description {
    recommended_name: Option<Name>,
    #[serde(default)]
    submission_names: Vec<Name>,
}
#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct Name {
    full_name: Value,
}
#[derive(Deserialize)]
struct Value {
    value: String,
}
#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct Gene {
    gene_name: Option<Value>,
    #[serde(default)]
    synonyms: Vec<Value>,
}
#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct Organism {
    scientific_name: String,
}
#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct Comment {
    comment_type: String,
    #[serde(default)]
    texts: Vec<Value>,
}
#[derive(Deserialize)]
struct Search {
    results: Vec<Record>,
}

pub fn normalize(record: Record, query: &str) -> ProteinResult {
    let function = record
        .comments
        .into_iter()
        .filter(|c| c.comment_type == "FUNCTION")
        .flat_map(|c| c.texts)
        .map(|t| t.value)
        .collect::<Vec<_>>()
        .join(" ");
    ProteinResult {
        query: query.into(),
        found: true,
        uniprot_id: Some(record.primary_accession),
        protein_name: record
            .protein_description
            .recommended_name
            .or_else(|| {
                record
                    .protein_description
                    .submission_names
                    .into_iter()
                    .next()
            })
            .map(|n| n.full_name.value),
        gene: record
            .genes
            .into_iter()
            .find_map(|g| g.gene_name.map(|v| v.value)),
        organism: Some(record.organism.scientific_name),
        function: (!function.is_empty()).then_some(function),
        go_annotations: vec![],
        cached: false,
        error: None,
    }
}
fn quote(value: &str) -> String {
    format!("\"{}\"", value.replace('\\', "\\\\").replace('"', "\\\""))
}
impl Upstream {
    pub async fn resolve(
        &self,
        query: &str,
        organism: &str,
    ) -> Result<Option<ProteinResult>, String> {
        for field in ["accession", "gene_exact"] {
            tracing::info!(event = "uniprot_lookup", strategy = field, query);
            let mut expression = format!("{field}:{}", quote(query));
            if !organism.is_empty() {
                expression.push_str(&format!(" AND organism_name:{}", quote(organism)));
            }
            let response = self
                .get(
                    &format!("{}/uniprotkb/search", self.uniprot),
                    &[
                        ("query", expression),
                        ("format", "json".into()),
                        ("size", "500".into()),
                    ],
                )
                .await?;
            if field == "accession" && matches!(response.status().as_u16(), 400 | 404) {
                continue;
            }
            if !response.status().is_success() {
                return Err(format!(
                    "UniProt returned HTTP {}",
                    response.status().as_u16()
                ));
            }
            // Never select from an incomplete set of candidates.
            if response
                .headers()
                .get("link")
                .and_then(|v| v.to_str().ok())
                .is_some_and(|v| v.contains("rel=\"next\""))
            {
                return Err("UniProt candidate set too large; specify an organism".into());
            }
            let search: Search = response
                .json()
                .await
                .map_err(|_| "invalid UniProt response")?;
            let mut exact: Vec<_> = search
                .results
                .into_iter()
                .filter(|r| {
                    (organism.is_empty()
                        || r.organism.scientific_name.eq_ignore_ascii_case(organism))
                        && if field == "accession" {
                            r.primary_accession.eq_ignore_ascii_case(query)
                        } else {
                            r.genes.iter().any(|g| {
                                g.gene_name
                                    .iter()
                                    .chain(g.synonyms.iter())
                                    .any(|n| n.value.eq_ignore_ascii_case(query))
                            })
                        }
                })
                .collect();
            if exact
                .iter()
                .any(|r| r.entry_type == "UniProtKB reviewed (Swiss-Prot)")
            {
                exact.retain(|r| r.entry_type == "UniProtKB reviewed (Swiss-Prot)");
            }
            match exact.len() {
                0 => continue,
                1 => return Ok(Some(normalize(exact.remove(0), query))),
                _ => {
                    return Err(
                        "ambiguous exact match; specify a more precise identifier or organism"
                            .into(),
                    );
                }
            }
        }
        Ok(None)
    }
}

/// Minimal organism metadata used by KEGG; independent of GO annotation lookup.
#[derive(Clone, Debug, serde::Serialize, Deserialize)]
pub struct ProteinTaxonomy {
    pub ncbi_taxonomy_id: Option<u64>,
    pub organism_name: Option<String>,
}
#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct TaxonomyRecord {
    primary_accession: String,
    organism: TaxonomyOrganism,
}
#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct TaxonomyOrganism {
    taxon_id: Option<u64>,
    scientific_name: Option<String>,
}
#[derive(Deserialize)]
struct TaxonomySearch {
    results: Vec<TaxonomyRecord>,
}
impl Upstream {
    /// Caller supplies validated canonical accessions, in bounded batches.
    pub async fn taxonomies(
        &self,
        accessions: &[String],
    ) -> Result<std::collections::BTreeMap<String, Option<ProteinTaxonomy>>, String> {
        let response = self
            .get(
                &format!("{}/uniprotkb/search", self.uniprot),
                &[
                    (
                        "query",
                        accessions
                            .iter()
                            .map(|id| format!("accession:{}", quote(id)))
                            .collect::<Vec<_>>()
                            .join(" OR "),
                    ),
                    ("fields", "accession,organism_id,organism_name".into()),
                    ("format", "json".into()),
                    ("size", "500".into()),
                ],
            )
            .await?;
        if !response.status().is_success() {
            return Err(format!(
                "UniProt taxonomy lookup returned HTTP {}",
                response.status().as_u16()
            ));
        }
        if response
            .headers()
            .get("link")
            .and_then(|v| v.to_str().ok())
            .is_some_and(|v| v.contains("rel=\"next\""))
        {
            return Err("UniProt taxonomy response is incomplete".into());
        }
        let search: TaxonomySearch = response
            .json()
            .await
            .map_err(|_| "invalid UniProt taxonomy response")?;
        let mut results: std::collections::BTreeMap<_, _> =
            accessions.iter().map(|id| (id.clone(), None)).collect();
        for record in search.results {
            // Accession queries may match secondary accessions. Do not assign
            // another primary accession's taxonomy to an unverified identifier.
            if let Some(value) = results.get_mut(&record.primary_accession) {
                if value.is_some() {
                    return Err("duplicate UniProt taxonomy record".into());
                }
                *value = Some(ProteinTaxonomy {
                    ncbi_taxonomy_id: record.organism.taxon_id.filter(|id| *id > 0),
                    organism_name: record.organism.scientific_name,
                });
            }
        }
        Ok(results)
    }
}
