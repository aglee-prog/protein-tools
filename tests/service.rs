use axum::{
    Json, Router,
    extract::{Query, State},
    http::StatusCode,
    routing::get,
};
use protein_tools::{
    api::App,
    cache::Cache,
    model::{ProteinRequest, ProteinResult, normalize},
    upstream::Upstream,
};
use serde_json::{Value, json};
use std::{
    collections::HashMap,
    sync::{Arc, Mutex},
};
use tokio::sync::Semaphore;

fn record(accession: &str, reviewed: bool) -> Value {
    json!({"primaryAccession":accession,"entryType":if reviewed {"UniProtKB reviewed (Swiss-Prot)"} else {"UniProtKB unreviewed (TrEMBL)"},"proteinDescription":{"recommendedName":{"fullName":{"value":"Cellular tumor antigen p53"}}},"genes":[{"geneName":{"value":"TP53"}}],"organism":{"scientificName":"Homo sapiens"},"comments":[{"commentType":"FUNCTION","texts":[{"value":"Authoritative function."}]}]})
}
type Calls = Arc<Mutex<Vec<String>>>;
async fn uniprot(
    State(calls): State<Calls>,
    Query(params): Query<HashMap<String, String>>,
) -> (StatusCode, Json<Value>) {
    let query = &params["query"];
    calls.lock().unwrap().push(query.clone());
    if query.contains("FAIL") {
        return (StatusCode::SERVICE_UNAVAILABLE, Json(json!({})));
    }
    if query.starts_with("accession:") {
        return (StatusCode::BAD_REQUEST, Json(json!({})));
    }
    let records = if query.contains("AMBIGUOUS") {
        let mut a = record("A", true);
        a["genes"][0]["geneName"]["value"] = json!("AMBIGUOUS");
        let mut b = a.clone();
        b["primaryAccession"] = json!("B");
        vec![a, b]
    } else if query.contains("MISSING") {
        vec![record("UNRELATED", true)]
    } else {
        vec![record("UNREVIEWED", false), record("P04637", true)]
    };
    (StatusCode::OK, Json(json!({"results": records})))
}
async fn quickgo(
    State(calls): State<Calls>,
    Query(params): Query<HashMap<String, String>>,
) -> Json<Value> {
    calls.lock().unwrap().push(format!(
        "quickgo:{}:{}",
        params["geneProductId"], params["page"]
    ));
    Json(
        json!({"pageInfo":{"current":params["page"].parse::<usize>().unwrap(),"total":2},"results":[{"geneProductId":"UniProtKB:P04637","goId":"GO:0003677","goAspect":"molecular_function","evidenceCode":"ECO:0000269","qualifier":"enables"}]}),
    )
}
async fn setup() -> (App, Calls, tempfile::TempDir, tokio::task::JoinHandle<()>) {
    let calls = Calls::default();
    let router = Router::new()
        .route("/uniprotkb/search", get(uniprot))
        .route("/annotation/search", get(quickgo))
        .with_state(calls.clone());
    let listener = tokio::net::TcpListener::bind("127.0.0.1:0").await.unwrap();
    let url = format!("http://{}", listener.local_addr().unwrap());
    let task = tokio::spawn(async move {
        axum::serve(listener, router).await.unwrap();
    });
    let dir = tempfile::tempdir().unwrap();
    let app = App {
        cache: Cache::open(dir.path().join("cache.sqlite").to_str().unwrap(), 30).unwrap(),
        upstream: Upstream::new(url.clone(), url).unwrap(),
        permits: Arc::new(Semaphore::new(4)),
    };
    (app, calls, dir, task)
}
#[test]
fn cache_key_normalization() {
    assert_eq!(
        normalize(" tp53 ", Some(" Homo sapiens ")),
        normalize("TP53", Some("homo SAPIENS"))
    );
    assert_eq!(normalize("EGFR", None), normalize("egfr", Some(" ")));
    assert_ne!(
        normalize("TP53", None),
        normalize("TP53", Some("Mus musculus"))
    );
}
#[tokio::test]
async fn cache_persistence_and_expiration() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("cache.sqlite");
    let path = path.to_str().unwrap();
    let key = normalize("TP53", None);
    let mut result = ProteinResult::missing("TP53", None);
    result.found = true;
    Cache::open(path, 30)
        .unwrap()
        .put(key.clone(), result.clone())
        .await
        .unwrap();
    assert!(
        Cache::open(path, 30)
            .unwrap()
            .get(key.clone())
            .await
            .unwrap()
            .is_some()
    );
    rusqlite::Connection::open(path)
        .unwrap()
        .execute("UPDATE cache SET saved=saved-2592001", [])
        .unwrap();
    assert!(
        Cache::open(path, 30)
            .unwrap()
            .get(key.clone())
            .await
            .unwrap()
            .is_none()
    );
    result.error = Some("network failure".into());
    Cache::open(path, 30)
        .unwrap()
        .put(key.clone(), result)
        .await
        .unwrap();
    assert!(
        Cache::open(path, 30)
            .unwrap()
            .get(key)
            .await
            .unwrap()
            .is_none()
    );
}
#[tokio::test]
async fn invalid_accession_falls_back_to_exact_gene_and_cache_skips_all_upstream() {
    let (app, calls, _dir, task) = setup().await;
    let result = app.lookup("TP53".into(), Some("Homo sapiens".into())).await;
    assert!(result.found);
    assert!(result.error.is_none());
    assert!(!result.cached);
    assert_eq!(result.uniprot_id.as_deref(), Some("P04637"));
    assert_eq!(result.gene.as_deref(), Some("TP53"));
    assert_eq!(result.organism.as_deref(), Some("Homo sapiens"));
    assert_eq!(
        result.protein_name.as_deref(),
        Some("Cellular tumor antigen p53")
    );
    assert_eq!(result.function.as_deref(), Some("Authoritative function."));
    assert_eq!(result.go_annotations.len(), 1);
    let before = calls.lock().unwrap().clone();
    assert!(before[0].starts_with("accession:"));
    assert!(before[1].starts_with("gene_exact:"));
    assert_eq!(before.len(), 4);
    let cached = app
        .lookup(" tp53 ".into(), Some(" homo SAPIENS ".into()))
        .await;
    assert!(cached.cached);
    assert_eq!(cached.query, " tp53 ");
    assert_eq!(*calls.lock().unwrap(), before);
    task.abort();
}
#[tokio::test]
async fn one_failure_does_not_invalidate_batch_and_unrelated_hits_are_rejected() {
    let (app, _calls, _dir, task) = setup().await;
    let results = app
        .batch(ProteinRequest {
            proteins: vec![
                "FAIL".into(),
                "TP53".into(),
                "MISSING".into(),
                "AMBIGUOUS".into(),
            ],
            organism: Some("Homo sapiens".into()),
        })
        .await;
    assert_eq!(results.len(), 4);
    assert!(!results[0].found);
    assert!(results[0].error.is_some());
    assert!(results[1].found);
    assert!(!results[2].found);
    assert!(results[2].error.is_none());
    assert!(!results[3].found);
    assert!(results[3].error.as_ref().unwrap().contains("ambiguous"));
    task.abort();
}
#[test]
fn go_deduplication_keeps_negation_and_filters_other_products() {
    let input = json!([
        {"geneProductId":"UniProtKB:P04637","goId":"GO:1","qualifier":"enables"},
        {"geneProductId":"UniProtKB:P04637","goId":"GO:1","qualifier":"enables"},
        {"geneProductId":"UniProtKB:P04637","goId":"GO:1","qualifier":"NOT|enables"},
        {"geneProductId":"UniProtKB:OTHER","goId":"GO:2"}]);
    let normalized =
        protein_tools::quickgo::normalize(serde_json::from_value(input).unwrap(), "P04637");
    assert_eq!(normalized.len(), 2);
}
