//! Compact cross-organism comparison over the existing KEGG lookup and cache.
use crate::{
    cache::Cache,
    kegg::{Kegg, KeggRequest, KeggStatus},
    upstream::Upstream,
};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use utoipa::ToSchema;

#[derive(Deserialize, ToSchema)]
pub struct KeggCompareRequest {
    /// Canonical UniProt accessions from any organisms; invalid entries produce per-input errors.
    #[schema(min_items = 2, max_items = 500)]
    pub proteins: Vec<String>,
    /// Minimum distinct accessions per pathway, default 2. Larger than the input count yields no matches.
    #[schema(minimum = 1, maximum = 500)]
    pub min_proteins: Option<usize>,
}
#[derive(Debug, Serialize, ToSchema)]
pub struct KeggCompareProtein {
    pub query: String,
    pub uniprot_id: String,
    pub status: KeggStatus,
    pub cached: bool,
    pub error: Option<String>,
}
#[derive(Debug, Serialize, ToSchema, PartialEq, Eq, PartialOrd, Ord)]
pub struct KeggCompareMember {
    pub uniprot_id: String,
    pub kegg_gene_id: String,
    pub pathway_id: String,
}
#[derive(Debug, Serialize, ToSchema)]
pub struct KeggComparePathway {
    /// Five-digit KEGG reference pathway number, preserving leading zeros.
    pub pathway_key: String,
    pub pathway_name: Option<String>,
    /// Unique accession/gene/pathway tuples, sorted lexically. Multiple genes may represent one protein.
    pub proteins: Vec<KeggCompareMember>,
    /// Distinct canonical accessions, not gene mappings or input occurrences.
    pub protein_count: usize,
}
#[derive(Debug, Serialize, ToSchema)]
pub struct KeggCompareResponse {
    pub complete: bool,
    /// Distinct valid canonical accessions supplied, including unmapped/failed accessions.
    pub protein_count: usize,
    pub matched_pathway_count: usize,
    /// Compact status per input, in input order, including duplicates and invalid inputs.
    pub proteins: Vec<KeggCompareProtein>,
    /// Only groups meeting the threshold, sorted by pathway key.
    pub pathways: Vec<KeggComparePathway>,
}

impl Kegg {
    pub async fn compare(
        &self,
        cache: &Cache,
        upstream: &Upstream,
        request: KeggCompareRequest,
    ) -> KeggCompareResponse {
        let minimum = request.min_proteins.unwrap_or(2);
        let lookup = self
            .lookup(
                cache,
                upstream,
                KeggRequest {
                    proteins: request.proteins,
                    pathway_id: None,
                    min_proteins: None,
                },
            )
            .await;
        let protein_count = lookup
            .proteins
            .iter()
            .filter_map(|p| crate::kegg::accession(&p.uniprot_id).ok())
            .collect::<BTreeSet<_>>()
            .len();
        let mut groups: BTreeMap<String, (BTreeSet<String>, BTreeSet<KeggCompareMember>)> =
            BTreeMap::new();
        for protein in &lookup.proteins {
            for path in &protein.pathways {
                // Use the resolved organism, never a fixed-width organism prefix or a name match.
                let Some(key) = protein.kegg_organism_codes.iter().find_map(|code| {
                    path.pathway_id
                        .strip_prefix(code)
                        .filter(|key| key.len() == 5 && key.bytes().all(|b| b.is_ascii_digit()))
                }) else {
                    continue;
                };
                let (names, members) = groups.entry(key.into()).or_default();
                if let Some(name) = &path.pathway_name {
                    // KEGG catalogs append organism metadata after the final " - ".
                    names.insert(
                        name.rsplit_once(" - ")
                            .map_or(name.as_str(), |(base, _)| base)
                            .trim()
                            .to_owned(),
                    );
                }
                members.insert(KeggCompareMember {
                    uniprot_id: protein.uniprot_id.clone(),
                    kegg_gene_id: path.kegg_gene_id.clone(),
                    pathway_id: path.pathway_id.clone(),
                });
            }
        }
        let pathways: Vec<_> = groups
            .into_iter()
            .filter_map(|(pathway_key, (names, members))| {
                let protein_count = members
                    .iter()
                    .map(|m| &m.uniprot_id)
                    .collect::<BTreeSet<_>>()
                    .len();
                (protein_count >= minimum).then(|| KeggComparePathway {
                    pathway_key,
                    pathway_name: names.into_iter().next(),
                    proteins: members.into_iter().collect(),
                    protein_count,
                })
            })
            .collect();
        KeggCompareResponse {
            complete: lookup.complete,
            protein_count,
            matched_pathway_count: pathways.len(),
            pathways,
            proteins: lookup
                .proteins
                .into_iter()
                .map(|p| KeggCompareProtein {
                    query: p.query,
                    uniprot_id: p.uniprot_id,
                    status: p.status,
                    cached: p.cached,
                    error: p.error,
                })
                .collect(),
        }
    }
}
