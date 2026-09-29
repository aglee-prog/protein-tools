use axum::{
    Router,
    extract::{Query, State},
    http::{HeaderMap, StatusCode},
    routing::get,
};
use protein_tools::{
    api::{App, router},
    cache::Cache,
    model::{ProteinResult, normalize},
    upstream::Upstream,
};
use serde_json::{Value, json};
use std::{
    collections::HashMap,
    sync::{Arc, Mutex},
};
use tokio::sync::Semaphore;

type Calls = Arc<Mutex<Vec<HashMap<String, String>>>>;
async fn serve(router: Router) -> (String, tokio::task::JoinHandle<()>) {
    let listener = tokio::net::TcpListener::bind("127.0.0.1:0").await.unwrap();
    let url = format!("http://{}", listener.local_addr().unwrap());
    let task = tokio::spawn(async move { axum::serve(listener, router).await.unwrap() });
    (url, task)
}
async fn background(
    State((calls, mode)): State<(Calls, String)>,
    Query(params): Query<HashMap<String, String>>,
    headers: HeaderMap,
) -> (StatusCode, HeaderMap, String) {
    calls.lock().unwrap().push(params.clone());
    assert_eq!(params["query"], "proteome:UP000005640 AND organism_id:9606");
    assert_eq!(params["fields"], "accession,go_p");
    let second = params.contains_key("cursor");
    let mut response_headers = HeaderMap::new();
    response_headers.insert("x-total-results", "20".parse().unwrap());
    response_headers.insert(
        "x-uniprot-release",
        if second && mode == "release" {
            "new"
        } else {
            "test"
        }
        .parse()
        .unwrap(),
    );
    if mode == "missing_total" {
        response_headers.remove("x-total-results");
    }
    if !second || mode == "loop" {
        let host = headers["host"].to_str().unwrap();
        response_headers.insert(
            "link",
            format!(
                "<http://{host}/uniprotkb/search?fields=accession,go_p&cursor=next>; rel=\"next\""
            )
            .parse()
            .unwrap(),
        );
    }
    if mode == "failure" && second {
        return (
            StatusCode::SERVICE_UNAVAILABLE,
            response_headers,
            "failure".into(),
        );
    }
    let mut body = "Entry\tGene Ontology (biological process)\n".to_owned();
    let range = if second && mode != "duplicate" {
        10..20
    } else {
        0..10
    };
    for i in range {
        if mode == "truncated" && i == 19 {
            continue;
        }
        let terms = if i < 7 && mode != "empty_terms" {
            (1..=12)
                .map(|j| format!("process {j} [GO:{j:07}]"))
                .collect::<Vec<_>>()
                .join("; ")
        } else {
            String::new()
        };
        if mode == "malformed" && i == 19 {
            body.push_str("bad row\n");
        } else {
            body.push_str(&format!("P{i:05}\t{terms}\n"));
        }
    }
    (StatusCode::OK, response_headers, body)
}
async fn setup(mode: &str) -> (App, Calls, tempfile::TempDir, tokio::task::JoinHandle<()>) {
    let calls: Calls = Arc::new(Mutex::new(Vec::new()));
    let (url, task) = serve(
        Router::new()
            .route("/uniprotkb/search", get(background))
            .with_state((calls.clone(), mode.to_owned())),
    )
    .await;
    let dir = tempfile::tempdir().unwrap();
    let cache = Cache::open(dir.path().join("cache.sqlite").to_str().unwrap(), 30).unwrap();
    for i in 0..21 {
        let id = format!("P{i:05}");
        let mut protein = ProteinResult::missing(&id, None);
        protein.found = true;
        protein.uniprot_id = Some(id.clone());
        protein.organism = Some("Homo sapiens".into());
        cache
            .put(normalize(&id, Some("Homo sapiens")), protein.clone())
            .await
            .unwrap();
        if i == 0 {
            cache
                .put(normalize("GENE0", Some("Homo sapiens")), protein)
                .await
                .unwrap();
        }
    }
    (
        App {
            cache,
            upstream: Upstream::new(url.clone(), url.clone()).unwrap(),
            kegg: protein_tools::kegg::Kegg::new(url).unwrap(),
            permits: Arc::new(Semaphore::new(4)),
        },
        calls,
        dir,
        task,
    )
}
#[tokio::test]
async fn api_statistics_compact_response_result_reuse_and_restart() {
    let (app, calls, dir, upstream_task) = setup("").await;
    let original = app.clone();
    let (url, task) = serve(router(app)).await;
    let client = reqwest::Client::new();
    let request =
        json!({"proteins":["P00000","P00001","P00002","P00003","P00004","GENE0"," p00000 "]});
    let response = client
        .post(format!("{url}/enrich-proteins"))
        .json(&request)
        .send()
        .await
        .unwrap();
    assert_eq!(response.status(), StatusCode::OK);
    let compact: Value = response.json().await.unwrap();
    assert_eq!(compact["input_count"], 5);
    assert_eq!(compact["background_count"], 20);
    assert_eq!(compact["tested_terms"], 12);
    assert_eq!(compact["significant_terms"], 12);
    assert_eq!(compact["top_terms"].as_array().unwrap().len(), 10);
    assert!(compact.get("items").is_none());
    assert!(compact["top_terms"][0].get("proteins").is_none());
    let expected = 21.0 / 15504.0;
    assert!((compact["top_terms"][0]["p_value"].as_f64().unwrap() - expected).abs() < 1e-12);
    let id = compact["result_id"].as_str().unwrap();
    let full: Value = client
        .get(format!("{url}/results/{id}"))
        .send()
        .await
        .unwrap()
        .json()
        .await
        .unwrap();
    assert_eq!(full["items"].as_array().unwrap().len(), 12);
    assert_eq!(full["items"][0]["proteins"].as_array().unwrap().len(), 5);
    assert_eq!(full["uniprot_release"], "test");
    assert_eq!(calls.lock().unwrap().len(), 2);
    task.abort();
    upstream_task.abort();
    let mut restarted = original;
    restarted.cache = Cache::open(dir.path().join("cache.sqlite").to_str().unwrap(), 30).unwrap();
    let (url, task) = serve(router(restarted.clone())).await;
    let response = client
        .post(format!("{url}/enrich-proteins"))
        .json(&json!({"result_id":id}))
        .send()
        .await
        .unwrap();
    assert_eq!(response.status(), StatusCode::OK);
    let repeated: Value = response.json().await.unwrap();
    assert_eq!(repeated["top_terms"], compact["top_terms"]);
    // Existing protein-info array is preserved; its header gives access to full cached data.
    let response = client
        .post(format!("{url}/protein-info"))
        .json(&json!({"proteins":["GENE0"],"organism":"Homo sapiens","include":["identity"]}))
        .send()
        .await
        .unwrap();
    let protein_id = response.headers()["x-result-id"]
        .to_str()
        .unwrap()
        .to_owned();
    assert!(response.json::<Value>().await.unwrap().is_array());
    assert_eq!(
        client
            .post(format!("{url}/enrich-proteins"))
            .json(&json!({"result_id":protein_id}))
            .send()
            .await
            .unwrap()
            .status(),
        StatusCode::OK
    );
    let doc: Value = client
        .get(format!("{url}/openapi.json"))
        .send()
        .await
        .unwrap()
        .json()
        .await
        .unwrap();
    assert_eq!(
        doc["paths"]["/enrich-proteins"]["post"]["operationId"],
        "enrich_proteins"
    );
    task.abort();
    let expired = Cache::open(dir.path().join("cache.sqlite").to_str().unwrap(), 0).unwrap();
    assert!(expired.get_analysis::<Value>(id).await.unwrap().is_none());
}
#[tokio::test]
async fn invalid_requests_and_missing_results() {
    let (app, calls, _dir, upstream_task) = setup("").await;
    let (url, task) = serve(router(app)).await;
    let client = reqwest::Client::new();
    for body in [
        json!({}),
        json!({"proteins":[]}),
        json!({"proteins":[" "]}),
        json!({"proteins":vec!["P00000";501]}),
        json!({"proteins":["P00000"],"result_id":"bad"}),
        json!({"proteins":["P00000"],"source":"reactome"}),
        json!({"proteins":["P00000"],"background":"input"}),
        json!({"result_id":"go_bp_human_proteome_v1"}),
    ] {
        assert_eq!(
            client
                .post(format!("{url}/enrich-proteins"))
                .json(&body)
                .send()
                .await
                .unwrap()
                .status(),
            StatusCode::BAD_REQUEST
        );
    }
    let id = "result_00000000000000000000000000000000";
    assert_eq!(
        client
            .post(format!("{url}/enrich-proteins"))
            .json(&json!({"result_id":id}))
            .send()
            .await
            .unwrap()
            .status(),
        StatusCode::NOT_FOUND
    );
    assert_eq!(
        client
            .get(format!("{url}/results/{id}"))
            .send()
            .await
            .unwrap()
            .status(),
        StatusCode::NOT_FOUND
    );
    assert!(calls.lock().unwrap().is_empty());
    assert_eq!(
        client
            .post(format!("{url}/enrich-proteins"))
            .json(&json!({"proteins":["P00020"]}))
            .send()
            .await
            .unwrap()
            .status(),
        StatusCode::BAD_REQUEST
    );
    task.abort();
    upstream_task.abort();
}
#[tokio::test]
async fn partial_backgrounds_never_produce_or_cache_results() {
    for mode in [
        "failure",
        "truncated",
        "release",
        "loop",
        "duplicate",
        "missing_total",
        "empty_terms",
        "malformed",
    ] {
        let (app, _calls, dir, upstream_task) = setup(mode).await;
        let (url, task) = serve(router(app)).await;
        let response = reqwest::Client::new()
            .post(format!("{url}/enrich-proteins"))
            .json(&json!({"proteins":["P00000"]}))
            .send()
            .await
            .unwrap();
        assert_eq!(response.status(), StatusCode::BAD_GATEWAY, "{mode}");
        let connection = rusqlite::Connection::open(dir.path().join("cache.sqlite")).unwrap();
        let count: i64 = connection
            .query_row("SELECT COUNT(*) FROM analysis_cache", [], |r| r.get(0))
            .unwrap();
        assert_eq!(count, 0, "{mode}");
        task.abort();
        upstream_task.abort();
    }
}

#[tokio::test]
async fn cached_nonhuman_and_unresolved_sets_are_not_reinterpreted() {
    let (app, calls, _dir, upstream_task) = setup("").await;
    let mut wrong_species = ProteinResult::missing("GENE0", None);
    wrong_species.found = true;
    wrong_species.organism = Some("Mus musculus".into());
    wrong_species.uniprot_id = Some("MOUSE".into());
    let missing = ProteinResult::missing("GENE0", None);
    let mut ids = Vec::new();
    for record in [wrong_species, missing] {
        ids.push(
            app.cache
                .put_analysis(None, json!({"proteins":["GENE0"], "items":[record]}))
                .await
                .unwrap(),
        );
    }
    let (url, task) = serve(router(app)).await;
    let client = reqwest::Client::new();
    for id in ids {
        let response = client
            .post(format!("{url}/enrich-proteins"))
            .json(&json!({"result_id":id}))
            .send()
            .await
            .unwrap();
        assert_eq!(response.status(), StatusCode::BAD_REQUEST);
        assert!(
            response
                .text()
                .await
                .unwrap()
                .contains("unresolved or nonhuman")
        );
    }
    assert!(calls.lock().unwrap().is_empty());
    task.abort();
    upstream_task.abort();
}
