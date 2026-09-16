use crate::{
    cache::Cache,
    model::{ProteinRequest, ProteinResult, normalize},
    upstream::Upstream,
};
use axum::{
    Json, Router,
    extract::State,
    http::StatusCode,
    routing::{get, post},
};
use futures::{StreamExt, stream};
use std::{collections::HashMap, sync::Arc};
use tokio::sync::Semaphore;
use utoipa::OpenApi;

pub const MAX_REQUEST_PROTEINS: usize = 500;
pub const MAX_CONCURRENT_LOOKUPS: usize = 4;

#[derive(Clone)]
pub struct App {
    pub cache: Cache,
    pub upstream: Upstream,
    pub permits: Arc<Semaphore>,
}
impl App {
    pub async fn lookup(&self, query: String, organism: Option<String>) -> ProteinResult {
        let key = normalize(&query, organism.as_deref());
        if key.0.is_empty() || key.0.len() > 128 || key.0.chars().any(char::is_control) {
            return ProteinResult::missing(&query, Some("invalid identifier".into()));
        }
        let _permit = match self.permits.acquire().await {
            Ok(p) => p,
            Err(_) => return ProteinResult::missing(&query, Some("server shutting down".into())),
        };
        match self.cache.get(key.clone()).await {
            Ok(Some(mut result)) => {
                tracing::debug!(event="cache_hit", query=%key.0);
                result.query = query;
                result.cached = true;
                return result;
            }
            Err(error) => tracing::warn!(event="cache_read_failure", %error),
            _ => {}
        }
        tracing::debug!(event="cache_miss", query=%key.0);
        let mut result = match self.upstream.resolve(&key.0, &key.1).await {
            Ok(Some(result)) => result,
            Ok(None) => return ProteinResult::missing(&query, None),
            Err(error) => return ProteinResult::missing(&query, Some(error)),
        };
        result.query = query;
        match self
            .upstream
            .annotations(result.uniprot_id.as_deref().expect("resolved accession"))
            .await
        {
            Ok(annotations) => result.go_annotations = annotations,
            Err(error) => {
                result.error = Some(error);
                return result;
            }
        }
        if let Err(error) = self.cache.put(key, result.clone()).await {
            tracing::warn!(event="cache_write_failure", %error);
        }
        result
    }
    pub async fn batch(&self, request: ProteinRequest) -> Vec<ProteinResult> {
        let mut indices = HashMap::new();
        let mut unique = Vec::new();
        let mut order = Vec::with_capacity(request.proteins.len());
        for query in &request.proteins {
            let key = normalize(query, request.organism.as_deref());
            let index = *indices.entry(key).or_insert_with(|| {
                unique.push(query.clone());
                unique.len() - 1
            });
            order.push(index);
        }
        let results: Vec<_> = stream::iter(
            unique
                .into_iter()
                .map(|query| self.lookup(query, request.organism.clone())),
        )
        .buffered(MAX_CONCURRENT_LOOKUPS)
        .collect()
        .await;
        let cache_hits = results.iter().filter(|result| result.cached).count();
        let resolved_count = order.iter().filter(|&&index| results[index].found).count();
        tracing::info!(
            event = "protein_request_summary",
            requested_count = order.len(),
            unique_count = results.len(),
            cache_hits,
            cache_misses = results.len() - cache_hits,
            resolved_count,
            unresolved_count = order.len() - resolved_count,
        );
        request
            .proteins
            .into_iter()
            .zip(order)
            .map(|(query, index)| {
                let mut result = results[index].clone();
                result.query = query;
                result
            })
            .collect()
    }
}
#[utoipa::path(post, path="/protein-info", operation_id="lookupProteinInfo", request_body=ProteinRequest,
    responses((status=200, description="One result per input, in input order. found indicates UniProt resolution; error may indicate incomplete QuickGO annotations. Only complete successful results are cached.", body=Vec<ProteinResult>), (status=400, description="Invalid batch size or organism")))]
/// Look up authoritative UniProt and Gene Ontology information for up to 500 proteins. Send the complete protein list in one request; clients should not manually split lists into multiple tool calls. The service handles batching, bounded concurrency, and caching internally. This service only retrieves facts; it does not interpret, classify, or group proteins. Only exact accessions or gene symbols are resolved; ambiguous matches are reported as unresolved.
async fn protein_info(
    State(app): State<App>,
    Json(request): Json<ProteinRequest>,
) -> Result<Json<Vec<ProteinResult>>, (StatusCode, &'static str)> {
    if request.proteins.is_empty() || request.proteins.len() > MAX_REQUEST_PROTEINS {
        return Err((
            StatusCode::BAD_REQUEST,
            "proteins must contain 1 to 500 identifiers",
        ));
    }
    if request
        .organism
        .as_ref()
        .is_some_and(|o| o.len() > 200 || o.chars().any(char::is_control))
    {
        return Err((StatusCode::BAD_REQUEST, "invalid organism"));
    }
    Ok(Json(app.batch(request).await))
}
#[utoipa::path(get, path="/health", operation_id="health", responses((status=200, description="Service is running; upstream reachability is not checked")))]
async fn health() -> Json<serde_json::Value> {
    Json(serde_json::json!({"status":"ok"}))
}
#[derive(OpenApi)]
#[openapi(
    paths(protein_info, health),
    components(schemas(ProteinRequest, ProteinResult, crate::model::GoAnnotation))
)]
pub struct ApiDoc;
pub fn router(app: App) -> Router {
    Router::new()
        .route("/protein-info", post(protein_info))
        .route("/health", get(health))
        .route("/openapi.json", get(|| async { Json(ApiDoc::openapi()) }))
        .with_state(app)
}
