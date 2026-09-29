use crate::enrichment::{EnrichmentRequest, EnrichmentResponse, valid_result_id};
use crate::enrichment::{
    EnrichmentTermDetail, EnrichmentTermPage, EnrichmentTermProteins, TermPageQuery, TermStatistics,
};
use crate::enrichment::{
    GroupEnrichmentDetail, GroupEnrichmentRequest, GroupEnrichmentResponse, GroupSelector,
};
use crate::kegg_compare::{KeggCompareRequest, KeggCompareResponse};
use crate::{
    cache::Cache,
    kegg::{Kegg, KeggRequest, KeggResponse, pathway_id},
    model::{ProteinRequest, ProteinResponse, ProteinResult, default_include, normalize},
    upstream::Upstream,
};
use axum::{
    Json, Router,
    extract::{Path, Query, State, rejection::JsonRejection},
    http::{HeaderMap, StatusCode},
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
    pub kegg: Kegg,
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
    responses((status=200, description="One result per input, in input order. query, found, cached, and error are always returned; selected fields are included when available. error may indicate incomplete QuickGO annotations.", body=Vec<ProteinResponse>, headers(("X-Result-Id" = String, description="Full cached result ID, when storage succeeds; accepted by enrich_proteins"))), (status=400, description="Invalid batch size, organism, include value, or GO limit")))]
/// Look up authoritative UniProt and Gene Ontology information for a protein set. Send the complete protein set in one request. Select only the information needed using `include` and control GO response size with `max_go_terms_per_protein`. The service handles batching, concurrency, and caching internally.
async fn protein_info(
    State(app): State<App>,
    request: Result<Json<ProteinRequest>, JsonRejection>,
) -> Result<(HeaderMap, Json<Vec<ProteinResponse>>), (StatusCode, String)> {
    let Json(request) = request.map_err(|error| (StatusCode::BAD_REQUEST, error.body_text()))?;
    if request.proteins.is_empty() || request.proteins.len() > MAX_REQUEST_PROTEINS {
        return Err((
            StatusCode::BAD_REQUEST,
            "proteins must contain 1 to 500 identifiers".into(),
        ));
    }
    if request
        .organism
        .as_ref()
        .is_some_and(|o| o.len() > 200 || o.chars().any(char::is_control))
    {
        return Err((StatusCode::BAD_REQUEST, "invalid organism".into()));
    }
    if request.max_go_terms_per_protein == Some(0) {
        return Err((
            StatusCode::BAD_REQUEST,
            "max_go_terms_per_protein must be at least 1".into(),
        ));
    }
    let include = request.include.clone().unwrap_or_else(default_include);
    let limit = request.max_go_terms_per_protein;
    let proteins = request.proteins.clone();
    let results = app.batch(request).await;
    let mut headers = HeaderMap::new();
    match app
        .cache
        .put_analysis(
            None,
            serde_json::json!({"proteins": proteins, "items": results}),
        )
        .await
    {
        Ok(id) => {
            headers.insert("x-result-id", id.parse().expect("generated result ID"));
        }
        Err(error) => tracing::warn!(event="result_cache_write_failure", %error),
    }
    Ok((
        headers,
        Json(
            results
                .into_iter()
                .map(|result| result.select(&include, limit))
                .collect(),
        ),
    ))
}
#[utoipa::path(post, path="/kegg-pathways", operation_id="lookupKeggPathways", request_body=KeggRequest,
    responses((status=200, description="Multi-organism KEGG associations and pathway groups. Check complete and per-protein errors before interpreting missing associations.", body=KeggResponse), (status=400, description="Invalid batch size, pathway ID, or minimum protein count")))]
/// Look up KEGG pathways for a mixed-organism set of UniProt accessions. Organisms are resolved automatically from UniProt taxonomy and KEGG. Returns per-protein mappings and pathway groups listing supplied members. Use min_proteins=2 for shared pathways or pathway_id to select membership in one pathway. Send the complete set; batching, rate limits, and persistent caching are automatic. Pathway groups retain their organism-specific IDs; no cross-species equivalence, enrichment, or relationships are inferred.
async fn kegg_pathways(
    State(app): State<App>,
    request: Result<Json<KeggRequest>, JsonRejection>,
) -> Result<Json<KeggResponse>, (StatusCode, String)> {
    let Json(mut request) =
        request.map_err(|error| (StatusCode::BAD_REQUEST, error.body_text()))?;
    if request.proteins.is_empty() || request.proteins.len() > MAX_REQUEST_PROTEINS {
        return Err((
            StatusCode::BAD_REQUEST,
            "proteins must contain 1 to 500 accessions".into(),
        ));
    }
    if request
        .min_proteins
        .is_some_and(|n| n == 0 || n > MAX_REQUEST_PROTEINS)
    {
        return Err((
            StatusCode::BAD_REQUEST,
            "min_proteins must be between 1 and 500".into(),
        ));
    }
    request.pathway_id = request
        .pathway_id
        .as_deref()
        .map(pathway_id)
        .transpose()
        .map_err(|error| (StatusCode::BAD_REQUEST, error))?;
    Ok(Json(
        app.kegg.lookup(&app.cache, &app.upstream, request).await,
    ))
}
#[utoipa::path(post, path="/compare-kegg-pathways", operation_id="compareKeggPathways", request_body=KeggCompareRequest,
    responses((status=200, description="Compact shared KEGG pathways across organisms, with per-input statuses. Check complete before interpreting missing matches.", body=KeggCompareResponse), (status=400, description="Invalid batch size or minimum protein count")))]
/// Compare KEGG pathways for 2–500 UniProt accessions across organisms. Groups by reference pathway number, defaults to at least 2 distinct proteins, and returns only matched pathways with provenance. Use min_proteins equal to the distinct input count for intersection. Reuses lookupKeggPathways provider and persistent cache.
async fn compare_kegg_pathways(
    State(app): State<App>,
    request: Result<Json<KeggCompareRequest>, JsonRejection>,
) -> Result<Json<KeggCompareResponse>, (StatusCode, String)> {
    let Json(request) = request.map_err(|error| (StatusCode::BAD_REQUEST, error.body_text()))?;
    if !(2..=MAX_REQUEST_PROTEINS).contains(&request.proteins.len()) {
        return Err((
            StatusCode::BAD_REQUEST,
            "proteins must contain 2 to 500 accessions".into(),
        ));
    }
    if request
        .min_proteins
        .is_some_and(|n| n == 0 || n > MAX_REQUEST_PROTEINS)
    {
        return Err((
            StatusCode::BAD_REQUEST,
            "min_proteins must be between 1 and 500".into(),
        ));
    }
    Ok(Json(
        app.kegg.compare(&app.cache, &app.upstream, request).await,
    ))
}
#[utoipa::path(post, path="/enrich-proteins", operation_id="enrich_proteins", request_body=EnrichmentRequest,
    responses((status=200, description="Compact GO BP overrepresentation result; complete data saved under result_id", body=EnrichmentResponse), (status=400, description="Invalid context, unresolved/ambiguous identifiers, or protein outside annotation universe"), (status=404, description="Result absent or expired"), (status=502, description="Upstream data unavailable or incomplete"), (status=500, description="Result cache or calculation failed")))]
/// Run deterministic GO BP enrichment with hypergeometric statistics and BH FDR. Supply proteins OR result_id and organism_taxon for new analyses. Returns a compact summary; full results stay behind result_id. Use get_enrichment_terms and get_enrichment_term for drill-down, then get_enrichment_term_proteins for selected hits. Do not retrieve the entire cached enrichment result. Cold background loading can take several minutes.
async fn enrich_proteins(
    State(app): State<App>,
    request: Result<Json<EnrichmentRequest>, JsonRejection>,
) -> Result<Json<EnrichmentResponse>, crate::enrichment::Error> {
    let Json(request) = request.map_err(|e| {
        (
            StatusCode::BAD_REQUEST,
            Json(serde_json::json!({"error":{"code":"invalid_request","message":e.body_text()}})),
        )
    })?;
    app.enrich(request).await.map(Json)
}
#[utoipa::path(post, path="/enrich-groups", operation_id="enrich_groups", request_body=GroupEnrichmentRequest,
    responses((status=200, description="Compact independent group summaries; check complete and failed_groups", body=GroupEnrichmentResponse), (status=400, description="Invalid request"), (status=502, description="Shared background unavailable"), (status=500, description="Cache or calculation failure")))]
/// Enrich 1–50 groups independently using one taxon/background and shared identifier resolution. At most 5,000 total inputs. Failed groups are explicit; valid groups remain available. No cross-group statistics. Inspect selected groups using get_group_enrichment and enrichment term tools with group selector.
async fn enrich_groups(
    State(app): State<App>,
    request: Result<Json<GroupEnrichmentRequest>, JsonRejection>,
) -> Result<Json<GroupEnrichmentResponse>, crate::enrichment::Error> {
    let Json(request) = request.map_err(|e| {
        (
            StatusCode::BAD_REQUEST,
            Json(serde_json::json!({"error":{"code":"invalid_request","message":e.body_text()}})),
        )
    })?;
    app.enrich_groups(request).await.map(Json)
}
#[utoipa::path(get, path="/results/{result_id}/enrichment/group", operation_id="get_group_enrichment", params(("result_id" = String, Path), ("group" = String, Query, description="Exact original group name")), responses((status=200, description="Selected group summary, provenance and explicit error if failed; no term collection", body=GroupEnrichmentDetail), (status=400, description="Invalid selector or result type"), (status=404, description="Result or group absent")))]
async fn get_group_enrichment(
    State(app): State<App>,
    Path(id): Path<String>,
    Query(query): Query<GroupSelector>,
) -> Result<Json<GroupEnrichmentDetail>, crate::enrichment::Error> {
    let group = query.group.ok_or_else(|| (StatusCode::BAD_REQUEST, Json(serde_json::json!({"error":{"code":"invalid_request","message":"group is required"}}))))?;
    app.group_enrichment_detail(id, group).await.map(Json)
}
#[utoipa::path(get, path="/results/{result_id}/enrichment/terms", operation_id="get_enrichment_terms", params(("result_id" = String, Path, description="Cached enrichment result ID"), ("group" = Option<String>, Query, description="Exact group name; required for batch results, omitted for single-set results"), ("limit" = Option<usize>, Query, description="Default 20; clamped to 50; minimum 1"), ("offset" = Option<usize>, Query, description="Default 0")),
    responses((status=200, description="Selective cached enrichment inspection", body=EnrichmentTermPage), (status=400, description="Invalid input or wrong result type"), (status=404, description="Result expired, absent, or term absent"), (status=422, description="Selected batch group failed; inspect get_group_enrichment for error details"), (status=500, description="Cache failure")))]
/// Browse significant terms in a bounded, paginated list. Ordered by FDR, p-value, then GO ID; includes nonsignificant terms. Limit defaults to 20 and is capped at 50.
async fn get_enrichment_terms(
    State(app): State<App>,
    Path(id): Path<String>,
    Query(query): Query<TermPageQuery>,
) -> Result<Json<EnrichmentTermPage>, crate::enrichment::Error> {
    let mut result = app
        .selected_enrichment_result(&id, query.group.as_deref())
        .await?;
    result.term_page(id, query).map(Json)
}
#[utoipa::path(get, path="/results/{result_id}/enrichment/terms/{term_id}", operation_id="get_enrichment_term", params(("result_id" = String, Path, description="Cached enrichment result ID"), ("group" = Option<String>, Query, description="Exact group name; required for batch results, omitted for single-set results"), ("term_id" = String, Path, description="GO term ID")),
    responses((status=200, description="Selective cached enrichment inspection", body=EnrichmentTermDetail), (status=400, description="Invalid input or wrong result type"), (status=404, description="Result expired, absent, or term absent"), (status=422, description="Selected batch group failed; inspect get_group_enrichment for error details"), (status=500, description="Cache failure")))]
/// Inspect statistics and small metadata for one selected GO term; no protein evidence or unrelated terms.
async fn get_enrichment_term(
    State(app): State<App>,
    Path((id, term_id)): Path<(String, String)>,
    Query(query): Query<GroupSelector>,
) -> Result<Json<EnrichmentTermDetail>, crate::enrichment::Error> {
    let result = app
        .selected_enrichment_result(&id, query.group.as_deref())
        .await?;
    let term = TermStatistics::from(result.selected_term(&term_id)?);
    Ok(Json(EnrichmentTermDetail {
        result_id: id,
        organism_taxon: result.organism_taxon,
        method: result.method,
        source: result.source,
        background: result.background,
        uniprot_release: result.uniprot_release,
        term,
    }))
}
#[utoipa::path(get, path="/results/{result_id}/enrichment/terms/{term_id}/proteins", operation_id="get_enrichment_term_proteins", params(("result_id" = String, Path, description="Cached enrichment result ID"), ("group" = Option<String>, Query, description="Exact group name; required for batch results, omitted for single-set results"), ("term_id" = String, Path, description="GO term ID")),
    responses((status=200, description="Selective cached enrichment inspection", body=EnrichmentTermProteins), (status=400, description="Invalid input or wrong result type"), (status=404, description="Result expired, absent, or term absent"), (status=422, description="Selected batch group failed; inspect get_group_enrichment for error details"), (status=500, description="Cache failure")))]
/// Use only after selecting a specific term. Returns its canonical input hits and original identifier mappings, with no additional UniProt or GO requests.
async fn get_enrichment_term_proteins(
    State(app): State<App>,
    Path((id, term_id)): Path<(String, String)>,
    Query(query): Query<GroupSelector>,
) -> Result<Json<EnrichmentTermProteins>, crate::enrichment::Error> {
    let result = app
        .selected_enrichment_result(&id, query.group.as_deref())
        .await?;
    let term = result.selected_term(&term_id)?;
    let identifier_mapping = result
        .identifier_mapping
        .iter()
        .filter(|m| term.proteins.contains(&m.canonical))
        .cloned()
        .collect();
    Ok(Json(EnrichmentTermProteins {
        result_id: id,
        term_id,
        hit_count: term.hit_count,
        proteins: term.proteins.clone(),
        identifier_mapping,
    }))
}
#[utoipa::path(get, path="/results/{result_id}", operation_id="getCachedResult", params(("result_id" = String, Path, description="Cached result ID")),
    responses((status=200, description="Complete cached result, or explicit truncated metadata for enrichment exceeding 32 KiB", body=serde_json::Value), (status=400, description="Invalid ID"), (status=404, description="Missing or expired result"), (status=500, description="Cache failure")))]
/// Retrieve cached data. Enrichment over 32 KiB returns truncated metadata only. Use selective enrichment tools for enrichment inspection.
async fn cached_result(
    State(app): State<App>,
    Path(id): Path<String>,
) -> Result<Json<serde_json::Value>, (StatusCode, String)> {
    if !valid_result_id(&id) {
        return Err((StatusCode::BAD_REQUEST, "invalid result_id".into()));
    }
    let value: serde_json::Value = app
        .cache
        .get_analysis(&id)
        .await
        .map_err(|e| (StatusCode::INTERNAL_SERVER_ERROR, e))?
        .ok_or((
            StatusCode::NOT_FOUND,
            "result_id not found or expired".into(),
        ))?;
    if crate::enrichment::is_enrichment(&value)
        && serde_json::to_vec(&value)
            .map_err(|e| (StatusCode::INTERNAL_SERVER_ERROR, e.to_string()))?
            .len()
            > crate::enrichment::MAX_DIRECT_ENRICHMENT_BYTES
    {
        if value["type"] == "group_enrichment" {
            return Ok(Json(serde_json::json!({
                "result_id": id, "type": "group_enrichment", "truncated": true,
                "group_count": value["groups"].as_object().map_or(0, |g| g.len()),
                "message": "Full result withheld: exceeds 32 KiB. Use get_group_enrichment or get_enrichment_terms / get_enrichment_term / get_enrichment_term_proteins with a group query selector."
            })));
        }
        return Ok(Json(
            serde_json::json!({"result_id": id, "type": "go_bp_enrichment",
            "item_count": value["items"].as_array().map_or(0, Vec::len), "truncated": true,
            "message": "Full result withheld: exceeds 32 KiB. Use get_enrichment_terms or get_enrichment_term, then get_enrichment_term_proteins for selected input hits."}),
        ));
    }
    Ok(Json(value))
}
#[utoipa::path(get, path="/health", operation_id="health", responses((status=200, description="Service is running; upstream reachability is not checked")))]
async fn health() -> Json<serde_json::Value> {
    Json(serde_json::json!({"status":"ok"}))
}
#[derive(OpenApi)]
#[openapi(
    paths(
        protein_info,
        kegg_pathways,
        compare_kegg_pathways,
        enrich_proteins,
        enrich_groups,
        get_group_enrichment,
        cached_result,
        get_enrichment_terms,
        get_enrichment_term,
        get_enrichment_term_proteins,
        health
    ),
    components(schemas(
        EnrichmentTermPage,
        EnrichmentTermDetail,
        EnrichmentTermProteins,
        GroupEnrichmentRequest,
        GroupEnrichmentResponse,
        GroupEnrichmentDetail,
        EnrichmentRequest,
        EnrichmentResponse,
        crate::enrichment::EnrichmentResult,
        crate::enrichment::EnrichmentTerm,
        crate::enrichment::TopTerm,
        KeggCompareRequest,
        KeggCompareResponse,
        crate::kegg_compare::KeggCompareProtein,
        crate::kegg_compare::KeggCompareMember,
        crate::kegg_compare::KeggComparePathway,
        KeggRequest,
        KeggResponse,
        crate::kegg::KeggProtein,
        crate::kegg::KeggAssociation,
        crate::kegg::KeggStatus,
        crate::kegg::KeggPathwayGroup,
        ProteinRequest,
        ProteinResponse,
        crate::model::Include,
        crate::model::GoAnnotation
    ))
)]
pub struct ApiDoc;
pub fn router(app: App) -> Router {
    Router::new()
        .route("/protein-info", post(protein_info))
        .route("/enrich-proteins", post(enrich_proteins))
        .route("/enrich-groups", post(enrich_groups))
        .route(
            "/results/{result_id}/enrichment/group",
            get(get_group_enrichment),
        )
        .route("/results/{result_id}", get(cached_result))
        .route(
            "/results/{result_id}/enrichment/terms",
            get(get_enrichment_terms),
        )
        .route(
            "/results/{result_id}/enrichment/terms/{term_id}",
            get(get_enrichment_term),
        )
        .route(
            "/results/{result_id}/enrichment/terms/{term_id}/proteins",
            get(get_enrichment_term_proteins),
        )
        .route("/kegg-pathways", post(kegg_pathways))
        .route("/compare-kegg-pathways", post(compare_kegg_pathways))
        .route("/health", get(health))
        .route("/openapi.json", get(|| async { Json(ApiDoc::openapi()) }))
        .with_state(app)
}
