use serde::{Deserialize, Serialize};
use utoipa::ToSchema;

#[derive(Deserialize, ToSchema)]
pub struct ProteinRequest {
    /// Gene symbols or UniProt accessions; send the complete set of 1 to 500 identifiers in one request.
    #[schema(min_items = 1, max_items = 500)]
    pub proteins: Vec<String>,
    /// Optional scientific organism name, for example Homo sapiens.
    pub organism: Option<String>,
}

#[derive(Clone, Debug, Deserialize, Serialize, ToSchema, PartialEq, Eq, PartialOrd, Ord)]
pub struct GoAnnotation {
    pub go_id: String,
    pub aspect: Option<String>,
    /// Evidence code supplied by QuickGO (usually an ECO identifier).
    pub evidence: Option<String>,
    /// Preserves negation and relation qualifiers from the authoritative annotation.
    pub qualifier: Option<String>,
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
