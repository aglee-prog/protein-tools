use crate::{model::GoAnnotation, upstream::Upstream};
use serde::Deserialize;
use std::collections::BTreeSet;
#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct Page {
    results: Vec<Annotation>,
    page_info: PageInfo,
}
#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct PageInfo {
    total: usize,
    current: usize,
}
#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Annotation {
    gene_product_id: String,
    go_id: String,
    go_aspect: Option<String>,
    evidence_code: Option<String>,
    qualifier: Option<String>,
}
pub fn normalize(annotations: Vec<Annotation>, accession: &str) -> Vec<GoAnnotation> {
    annotations
        .into_iter()
        .filter(|a| {
            a.gene_product_id == format!("UniProtKB:{accession}") || a.gene_product_id == accession
        })
        .map(|a| GoAnnotation {
            go_id: a.go_id,
            aspect: a.go_aspect,
            evidence: a.evidence_code,
            qualifier: a.qualifier,
        })
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect()
}
impl Upstream {
    pub async fn annotations(&self, accession: &str) -> Result<Vec<GoAnnotation>, String> {
        let mut annotations = BTreeSet::new();
        for page in 1..=100 {
            tracing::info!(event = "quickgo_lookup", accession, page);
            let response = self
                .get(
                    &format!("{}/annotation/search", self.quickgo),
                    &[
                        ("geneProductId", format!("UniProtKB:{accession}")),
                        ("limit", "200".into()),
                        ("page", page.to_string()),
                    ],
                )
                .await?;
            if !response.status().is_success() {
                return Err(format!(
                    "QuickGO returned HTTP {}",
                    response.status().as_u16()
                ));
            }
            let data: Page = response
                .json()
                .await
                .map_err(|_| "invalid QuickGO response")?;
            if data.page_info.current != page {
                return Err("unexpected QuickGO pagination".into());
            }
            annotations.extend(normalize(data.results, accession));
            if page >= data.page_info.total {
                return Ok(annotations.into_iter().collect());
            }
        }
        Err("QuickGO pagination limit exceeded; annotations incomplete".into())
    }
}
