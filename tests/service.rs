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
        let identifier = query.split('"').nth(1).unwrap();
        if identifier.starts_with('P') && identifier[1..].chars().all(|c| c.is_ascii_digit()) {
            return (
                StatusCode::OK,
                Json(json!({"results": [record(identifier, true)]})),
            );
        }
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

async fn serve_app(app: App) -> (String, tokio::task::JoinHandle<()>) {
    let listener = tokio::net::TcpListener::bind("127.0.0.1:0").await.unwrap();
    let url = format!("http://{}", listener.local_addr().unwrap());
    let task = tokio::spawn(async move {
        axum::serve(listener, protein_tools::api::router(app))
            .await
            .unwrap();
    });
    (url, task)
}

#[tokio::test]
async fn public_limits_large_mixed_batches_and_cache() {
    let (app, calls, _dir, upstream_task) = setup().await;
    let (url, task) = serve_app(app).await;
    let client = reqwest::Client::new();
    for count in [0, 501] {
        let response = client
            .post(format!("{url}/protein-info"))
            .json(&json!({"proteins":vec!["TP53"; count]}))
            .send()
            .await
            .unwrap();
        assert_eq!(response.status(), StatusCode::BAD_REQUEST);
        assert!(response.text().await.unwrap().contains("1 to 500"));
    }
    assert!(calls.lock().unwrap().is_empty());
    for count in [70, 200, 230, 500] {
        let proteins: Vec<_> = (0..count).map(|i| format!("P{i:05}")).collect();
        let before = calls.lock().unwrap().len();
        let response = client
            .post(format!("{url}/protein-info"))
            .json(&json!({"proteins": proteins, "organism":"Homo sapiens"}))
            .send()
            .await
            .unwrap();
        assert_eq!(response.status(), StatusCode::OK);
        let results: Vec<ProteinResult> = response.json().await.unwrap();
        assert_eq!(
            results.iter().map(|r| &r.query).collect::<Vec<_>>(),
            proteins.iter().collect::<Vec<_>>()
        );
        assert!(results.iter().all(|r| r.found && r.error.is_none()));
        let cached = match count {
            70 => 0,
            200 => 70,
            230 => 200,
            500 => 230,
            _ => unreachable!(),
        };
        assert_eq!(results.iter().filter(|r| r.cached).count(), cached);
        // One accession lookup and two QuickGO pages per uncached identifier.
        assert_eq!(calls.lock().unwrap().len() - before, (count - cached) * 3);
        let before = calls.lock().unwrap().len();
        let results: Vec<ProteinResult> = client
            .post(format!("{url}/protein-info"))
            .json(&json!({"proteins": proteins, "organism":"Homo sapiens"}))
            .send()
            .await
            .unwrap()
            .json()
            .await
            .unwrap();
        assert!(results.iter().all(|r| r.cached));
        assert_eq!(calls.lock().unwrap().len(), before);
    }
    let schema: Value = client
        .get(format!("{url}/openapi.json"))
        .send()
        .await
        .unwrap()
        .json()
        .await
        .unwrap();
    assert_eq!(
        schema["components"]["schemas"]["ProteinRequest"]["properties"]["proteins"]["maxItems"],
        500
    );
    assert!(
        schema["paths"]["/protein-info"]["post"]
            .to_string()
            .contains("complete protein list")
    );
    task.abort();
    upstream_task.abort();
}

#[tokio::test]
async fn duplicates_share_successes_and_failures_and_keep_original_queries() {
    let (app, calls, _dir, task) = setup().await;
    let proteins = vec!["TP53", "FAIL", " tp53 ", "MISSING", "fail", "missing"];
    let results = app
        .batch(ProteinRequest {
            proteins: proteins.iter().map(|s| (*s).into()).collect(),
            organism: None,
        })
        .await;
    assert_eq!(
        results.iter().map(|r| r.query.as_str()).collect::<Vec<_>>(),
        proteins
    );
    assert!(results[0].found && results[2].found);
    assert!(results[1].error.is_some() && results[4].error.is_some());
    assert!(!results[3].found && !results[5].found);
    assert!(results.iter().all(|r| !r.cached));
    assert_eq!(calls.lock().unwrap().len(), 7);
    task.abort();
}

#[tokio::test]
async fn concurrent_completion_keeps_order_and_shared_limit() {
    use std::sync::atomic::{AtomicUsize, Ordering};
    #[derive(Clone, Default)]
    struct Mock {
        active: Arc<AtomicUsize>,
        peak: Arc<AtomicUsize>,
        finished: Arc<Mutex<Vec<String>>>,
        gate: Arc<tokio::sync::Notify>,
    }
    let mock = Mock::default();
    let state = mock.clone();
    let router = Router::new().route("/uniprotkb/search", get(
        |State(mock): State<Mock>, Query(params): Query<HashMap<String, String>>| async move {
            let query = params["query"].split('"').nth(1).unwrap().to_string();
            let active = mock.active.fetch_add(1, Ordering::SeqCst) + 1;
            mock.peak.fetch_max(active, Ordering::SeqCst);
            if query == "P00000" {
                mock.gate.notified().await;
            } else {
                tokio::time::sleep(std::time::Duration::from_millis(5)).await;
            }
            mock.finished.lock().unwrap().push(query.clone());
            if query == "P00001" { mock.gate.notify_one(); }
            mock.active.fetch_sub(1, Ordering::SeqCst);
            // Isolated upstream failures also preserve the response shape.
            (StatusCode::SERVICE_UNAVAILABLE, Json(json!({})))
        }
    )).with_state(state);
    let listener = tokio::net::TcpListener::bind("127.0.0.1:0").await.unwrap();
    let url = format!("http://{}", listener.local_addr().unwrap());
    let task = tokio::spawn(async move {
        axum::serve(listener, router).await.unwrap();
    });
    let (mut app, _calls, _dir, other_task) = setup().await;
    app.upstream = Upstream::new(url.clone(), url).unwrap();
    let request = |start| ProteinRequest {
        proteins: (start..start + 70).map(|i| format!("P{i:05}")).collect(),
        organism: None,
    };
    let (a, b) = tokio::time::timeout(std::time::Duration::from_secs(10), async {
        tokio::join!(app.batch(request(0)), app.batch(request(70)))
    })
    .await
    .unwrap();
    for (results, start) in [(a, 0), (b, 70)] {
        assert_eq!(
            results.iter().map(|r| r.query.clone()).collect::<Vec<_>>(),
            request(start).proteins
        );
        assert!(results.iter().all(|r| !r.found && r.error.is_some()));
    }
    let finished = mock.finished.lock().unwrap();
    assert!(
        finished.iter().position(|q| q == "P00001").unwrap()
            < finished.iter().position(|q| q == "P00000").unwrap()
    );
    assert!((2..=4).contains(&mock.peak.load(Ordering::SeqCst)));
    task.abort();
    other_task.abort();
}
