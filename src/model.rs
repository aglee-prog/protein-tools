use serde::{Deserialize, Serialize};
use utoipa::ToSchema;

/// Allowed response categories: identity provides UniProt identifiers, name, gene, and organism; function provides the UniProt function comment; go.biological_process, go.molecular_function, and go.cellular_component select annotations for those GO aspects; go.evidence adds available QuickGO evidence codes to the selected GO annotations.
#[derive(Clone, Copy, Debug, Deserialize, ToSchema, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum Include {
    /// Identifiers, name, gene, and organism from UniProt.
    Identity,
    /// Normalized UniProt function comment.
    Function,
    #[serde(rename = "go.biological_process")]
    /// Biological Process GO annotations.
    GoBiologicalProcess,
    #[serde(rename = "go.molecular_function")]
    /// Molecular Function GO annotations.
    GoMolecularFunction,
    #[serde(rename = "go.cellular_component")]
    /// Cellular Component GO annotations.
    GoCellularComponent,
    #[serde(rename = "go.evidence")]
    /// Add available QuickGO evidence codes to selected GO annotations.
    GoEvidence,
}

pub fn default_include() -> Vec<Include> {
    vec![
        Include::Identity,
        Include::Function,
        Include::GoBiologicalProcess,
    ]
}

#[derive(Deserialize, ToSchema)]
pub struct ProteinRequest {
    /// Gene symbols or UniProt accessions; send the complete set of 1 to 500 identifiers in one request.
    #[schema(min_items = 1, max_items = 500)]
    pub proteins: Vec<String>,
    /// Optional scientific organism name, for example Homo sapiens.
    pub organism: Option<String>,
    /// Categories to return. Defaults to identity, function, and Biological Process. go.evidence augments selected GO aspects.
    pub include: Option<Vec<Include>>,
    /// Optional maximum total GO annotations per protein across selected aspects. Omit for no GO term limit. Explicit values must be at least 1; annotations are sorted deterministically before truncation.
    #[schema(minimum = 1)]
    pub max_go_terms_per_protein: Option<usize>,
}

#[derive(Clone, Debug, Deserialize, Serialize, ToSchema, PartialEq, Eq, PartialOrd, Ord)]
pub struct GoAnnotation {
    pub go_id: String,
    pub aspect: Option<String>,
    /// Evidence code supplied by QuickGO (usually an ECO identifier).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub evidence: Option<String>,
    /// Preserves negation and relation qualifiers from the authoritative annotation.
    pub qualifier: Option<String>,
}

#[derive(Clone, Debug, Serialize, ToSchema)]
pub struct ProteinResponse {
    pub query: String,
    pub found: bool,
    pub cached: bool,
    pub error: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub uniprot_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub protein_name: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gene: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub organism: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub function: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub go_annotations: Option<Vec<GoAnnotation>>,
}

impl ProteinResult {
    pub fn select(self, include: &[Include], limit: Option<usize>) -> ProteinResponse {
        let identity = include.contains(&Include::Identity);
        let aspects: Vec<&str> = include
            .iter()
            .filter_map(|field| match field {
                Include::GoBiologicalProcess => Some("biological_process"),
                Include::GoMolecularFunction => Some("molecular_function"),
                Include::GoCellularComponent => Some("cellular_component"),
                _ => None,
            })
            .collect();
        let go_annotations = (!aspects.is_empty()).then(|| {
            let mut annotations: Vec<_> = self
                .go_annotations
                .into_iter()
                .filter(|a| {
                    a.aspect
                        .as_deref()
                        .is_some_and(|aspect| aspects.contains(&aspect))
                })
                .collect();
            // Stable presentation order; no biological priority is assigned.
            annotations.sort_by(|a, b| {
                (&a.aspect, &a.go_id, &a.qualifier, &a.evidence).cmp(&(
                    &b.aspect,
                    &b.go_id,
                    &b.qualifier,
                    &b.evidence,
                ))
            });
            if let Some(limit) = limit {
                annotations.truncate(limit);
            }
            if !include.contains(&Include::GoEvidence) {
                for annotation in &mut annotations {
                    annotation.evidence = None;
                }
            }
            annotations
        });
        ProteinResponse {
            query: self.query,
            found: self.found,
            cached: self.cached,
            error: self.error,
            uniprot_id: identity.then_some(self.uniprot_id).flatten(),
            protein_name: identity.then_some(self.protein_name).flatten(),
            gene: identity.then_some(self.gene).flatten(),
            organism: identity.then_some(self.organism).flatten(),
            function: include
                .contains(&Include::Function)
                .then_some(self.function)
                .flatten(),
            go_annotations,
        }
    }
}

#[derive(Clone, Debug, Deserialize, Serialize, ToSchema)]
pub struct ProteinResult {
    pub query: String,
    pub found: bool,
    pub uniprot_id: Option<String>,
    pub protein_name: Option<String>,
    pub gene: Option<String>,
    pub organism: Option<String>,
    pub function: Option<String>,
    pub go_annotations: Vec<GoAnnotation>,
    pub cached: bool,
    pub error: Option<String>,
}
impl ProteinResult {
    pub fn missing(query: &str, error: Option<String>) -> Self {
        Self {
            query: query.into(),
            found: false,
            uniprot_id: None,
            protein_name: None,
            gene: None,
            organism: None,
            function: None,
            go_annotations: vec![],
            cached: false,
            error,
        }
    }
}
pub fn normalize(query: &str, organism: Option<&str>) -> (String, String) {
    (
        query.trim().to_ascii_uppercase(),
        organism.unwrap_or("").trim().to_lowercase(),
    )
}
