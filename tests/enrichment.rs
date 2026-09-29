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
    assert_eq!(params["fields"], "accession,reviewed,xref_geneid,go_p");
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
    let mut body = "Entry\tReviewed\tGeneID\tGene Ontology (biological process)\n".to_owned();
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
        } else if i < 19 && mode != "empty_terms" {
            "other [GO:0099999]".into()
        } else {
            String::new()
        };
        if mode == "malformed" && i == 19 {
            body.push_str("bad row\n");
        } else {
            body.push_str(&format!("P{i:05}\treviewed\t{i};\t{terms}\n"));
        }
    }
    (StatusCode::OK, response_headers, body)
}
async fn setup(mode: &str) -> (App, Calls, tempfile::TempDir, tokio::task::JoinHandle<()>) {
    let calls: Calls = Arc::new(Mutex::new(Vec::new()));
    let (url, task) = serve(
        Router::new()
            .route("/proteomes/search", get(proteome_fixture))
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
            .put(
                (id.clone(), "taxon:9606:resolver_v2".into()),
                protein.clone(),
            )
            .await
            .unwrap();
        if i == 0 {
            cache
                .put(
                    ("GENE0".into(), "taxon:9606:resolver_v2".into()),
                    protein.clone(),
                )
                .await
                .unwrap();
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
    let request = json!({"organism_taxon":9606,"proteins":["P00000","P00001","P00002","P00003","P00004","GENE0"," p00000 "]});
    let response = client
        .post(format!("{url}/enrich-proteins"))
        .json(&request)
        .send()
        .await
        .unwrap();
    assert_eq!(response.status(), StatusCode::OK);
    let compact: Value = response.json().await.unwrap();
    assert_eq!(compact["input_count"], 5);
    assert_eq!(compact["background_count"], 19);
    assert_eq!(compact["tested_terms"], 13);
    assert_eq!(compact["significant_terms"], 12);
    assert_eq!(compact["top_terms"].as_array().unwrap().len(), 10);
    assert!(compact.get("items").is_none());
    assert!(compact["top_terms"][0].get("proteins").is_none());
    let expected = 21.0 / 11628.0;
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
    assert_eq!(full["items"].as_array().unwrap().len(), 13);
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
        .json(&json!({"organism_taxon":9606,"proteins":["GENE0"],"organism":"Homo sapiens","include":["identity"]}))
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
            .json(&json!({"result_id":protein_id,"organism_taxon":9606}))
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
        json!({"organism_taxon":9606,"proteins":[]}),
        json!({"organism_taxon":9606,"proteins":[" "]}),
        json!({"organism_taxon":9606,"proteins":vec!["P00000";501]}),
        json!({"organism_taxon":9606,"proteins":["P00000"],"result_id":"bad"}),
        json!({"organism_taxon":9606,"proteins":["P00000"],"source":"reactome"}),
        json!({"organism_taxon":9606,"proteins":["P00000"],"background":"input"}),
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
            .json(&json!({"organism_taxon":9606,"proteins":["P00020"]}))
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
            .json(&json!({"organism_taxon":9606,"proteins":["P00000"]}))
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
async fn missing_context_and_cached_unresolved_are_structured() {
    let (app, calls, _dir, upstream_task) = setup("").await;
    let id = app
        .cache
        .put_analysis(None, json!({"items":[ProteinResult::missing("ATR", None)]}))
        .await
        .unwrap();
    let (url, task) = serve(router(app)).await;
    let client = reqwest::Client::new();
    for (request, code) in [
        (json!({"proteins":["ATR"]}), "organism_required"),
        (
            json!({"result_id":id,"organism_taxon":9606}),
            "resolution_failed",
        ),
        (
            json!({"proteins":["ATR"],"organism_taxon":0}),
            "organism_required",
        ),
    ] {
        let response = client
            .post(format!("{url}/enrich-proteins"))
            .json(&request)
            .send()
            .await
            .unwrap();
        assert_eq!(response.status(), StatusCode::BAD_REQUEST);
        assert_eq!(
            response.json::<Value>().await.unwrap()["error"]["code"],
            code
        );
    }
    assert!(calls.lock().unwrap().is_empty());
    task.abort();
    upstream_task.abort();
}

const LIVE_SYMBOLS: [&str; 10] = [
    "TP53", "ATM", "ATR", "BRCA1", "BRCA2", "CHEK1", "CHEK2", "RAD51", "MDM2", "CDKN1A",
];
const LIVE_ACCESSIONS: [&str; 10] = [
    "P04637", "Q13315", "Q13535", "P38398", "P51587", "O14757", "O96017", "Q06609", "Q00987",
    "P38936",
];
fn candidate(symbol: &str, accession: &str, taxon: u64) -> Value {
    json!({"primaryAccession":accession,"entryType":"UniProtKB reviewed (Swiss-Prot)",
        "proteinDescription":{},"genes":[{"geneName":{"value":symbol}}],
        "organism":{"scientificName":if taxon == 9606 {"Homo sapiens"} else {"Mus musculus"}, "taxonId":taxon}})
}
async fn resolution_fixture(
    State(calls): State<Calls>,
    Query(params): Query<HashMap<String, String>>,
) -> (StatusCode, HeaderMap, String) {
    calls.lock().unwrap().push(params.clone());
    let query = &params["query"];
    let taxon = if query.ends_with("organism_id:10090") {
        10090
    } else {
        9606
    };
    assert!(query.ends_with(&format!("organism_id:{taxon}")));
    if params["format"] == "tsv" {
        assert_eq!(
            query,
            &format!("proteome:UP000005640 AND organism_id:{taxon}")
        );
        let mut headers = HeaderMap::new();
        headers.insert("x-total-results", "12".parse().unwrap());
        headers.insert("x-uniprot-release", "fixture".parse().unwrap());
        let mut body = "Entry\tReviewed\tGeneID\tGene Ontology (biological process)\n".to_owned();
        for (i, id) in LIVE_ACCESSIONS.iter().enumerate() {
            body.push_str(&format!(
                "{id}\treviewed\t{i};\tDNA damage response [GO:0006974]\n"
            ));
        }
        body.push_str(
            "Q00001\tunreviewed\t0;\tother [GO:0000001]\nQ00002\tunreviewed\t\troot [GO:0008150]\n",
        );
        return (StatusCode::OK, headers, body);
    }
    let identifier = query.split('"').nth(1).unwrap();
    let mut records = vec![];
    if query.starts_with("accession:") {
        if let Some(i) = LIVE_ACCESSIONS.iter().position(|id| *id == identifier) {
            records.push(candidate(LIVE_SYMBOLS[i], identifier, taxon));
        } else if ["Q00001", "Q00002"].contains(&identifier) {
            records.push(candidate("EXTRA", identifier, taxon));
        }
    } else if identifier == "AMBIGUOUS" {
        records.extend([
            candidate(identifier, "Q11111", taxon),
            candidate(identifier, "Q22222", taxon),
        ]);
    } else if identifier == "ALIAS" {
        let mut record = candidate("PRIMARY", "Q13535", taxon);
        record["genes"][0]["synonyms"] = json!([{"value":"ALIAS"}]);
        records.push(record);
    } else if identifier == "WRONGTAXON" {
        records.push(candidate(identifier, "Q13535", 10090));
    } else if let Some(i) = LIVE_SYMBOLS.iter().position(|id| *id == identifier) {
        records.push(candidate(identifier, LIVE_ACCESSIONS[i], taxon));
        // Actual ATR alias collision: both of these are reviewed human entries.
        if identifier == "ATR" {
            for (symbol, id) in [("SERPINA2", "P20848"), ("ANTXR1", "Q9H6X2")] {
                let mut record = candidate(symbol, id, taxon);
                record["genes"][0]["synonyms"] = json!([{"value":"ATR"}]);
                records.push(record);
            }
        }
        let mut unreviewed = candidate(identifier, "A0A0000001", taxon);
        unreviewed["entryType"] = json!("UniProtKB unreviewed (TrEMBL)");
        records.push(unreviewed);
    }
    (
        StatusCode::OK,
        HeaderMap::new(),
        json!({"results":records}).to_string(),
    )
}
async fn resolution_setup() -> (App, Calls, tempfile::TempDir, tokio::task::JoinHandle<()>) {
    let calls = Calls::default();
    let (url, task) = serve(
        Router::new()
            .route("/proteomes/search", get(proteome_fixture))
            .route("/uniprotkb/search", get(resolution_fixture))
            .with_state(calls.clone()),
    )
    .await;
    let dir = tempfile::tempdir().unwrap();
    (
        App {
            cache: Cache::open(dir.path().join("test.sqlite").to_str().unwrap(), 30).unwrap(),
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
async fn exact_live_set_resolves_all_ten_and_preserves_provenance() {
    let (app, calls, _dir, upstream_task) = resolution_setup().await;
    let (url, task) = serve(router(app)).await;
    let client = reqwest::Client::new();
    let response = client.post(format!("{url}/enrich-proteins"))
        .json(&json!({"proteins":LIVE_SYMBOLS,"organism_taxon":9606,"source":"go_bp","background":"proteome"}))
        .send().await.unwrap();
    assert_eq!(response.status(), StatusCode::OK);
    let compact: Value = response.json().await.unwrap();
    assert_eq!(compact["input_count"], 10);
    assert_eq!(compact["background_count"], 11);
    assert_eq!(compact["organism_taxon"], 9606);
    let full: Value = client
        .get(format!(
            "{url}/results/{}",
            compact["result_id"].as_str().unwrap()
        ))
        .send()
        .await
        .unwrap()
        .json()
        .await
        .unwrap();
    assert_eq!(full["input_identifiers"], json!(LIVE_SYMBOLS));
    for (i, id) in LIVE_ACCESSIONS.iter().enumerate() {
        assert_eq!(
            full["identifier_mapping"][i],
            json!({"input":LIVE_SYMBOLS[i],"canonical":id})
        );
        assert!(full["proteins"].as_array().unwrap().contains(&json!(id)));
    }
    assert_eq!(full["universe"]["proteome_entry_count"], 12);
    assert_eq!(full["universe"]["annotated_entry_count"], 11);
    assert_eq!(full["universe"]["unannotated_entry_count"], 1);
    assert_eq!(full["universe"]["unique_gene_id_count"], 10);
    assert_eq!(full["universe"]["entries_without_gene_id"], 1);
    assert_eq!(
        full["universe"]["eligible_accessions"]
            .as_array()
            .unwrap()
            .len(),
        11
    );
    assert!(
        calls
            .lock()
            .unwrap()
            .iter()
            .any(|p| p["query"] == "gene_exact:\"ATR\" AND organism_id:9606")
    );
    let response = client
        .post(format!("{url}/enrich-proteins"))
        .json(&json!({"result_id":compact["result_id"],"organism_taxon":10090}))
        .send()
        .await
        .unwrap();
    assert_eq!(response.status(), StatusCode::BAD_REQUEST);
    let response: Value = client
        .post(format!("{url}/enrich-proteins"))
        .json(&json!({"result_id":compact["result_id"]}))
        .send()
        .await
        .unwrap()
        .json()
        .await
        .unwrap();
    let repeated: Value = client
        .get(format!(
            "{url}/results/{}",
            response["result_id"].as_str().unwrap()
        ))
        .send()
        .await
        .unwrap()
        .json()
        .await
        .unwrap();
    assert_eq!(repeated["identifier_mapping"], full["identifier_mapping"]);
    task.abort();
    upstream_task.abort();
}
#[tokio::test]
async fn failures_are_aggregated_and_taxon_universes_are_isolated() {
    let (app, calls, _dir, upstream_task) = resolution_setup().await;
    let (url, task) = serve(router(app.clone())).await;
    let client = reqwest::Client::new();
    let response = client
        .post(format!("{url}/enrich-proteins"))
        .json(
            &json!({"proteins":["TP53","MISSING","AMBIGUOUS","WRONGTAXON"],"organism_taxon":9606}),
        )
        .send()
        .await
        .unwrap();
    assert_eq!(response.status(), StatusCode::BAD_REQUEST);
    let error: Value = response.json().await.unwrap();
    assert_eq!(error["error"]["code"], "resolution_failed");
    let issues = error["error"]["identifiers"].as_array().unwrap();
    assert_eq!(issues.len(), 3);
    assert!(
        issues
            .iter()
            .any(|p| p["identifier"] == "AMBIGUOUS" && p["reason"] == "ambiguous")
    );
    assert!(!calls.lock().unwrap().iter().any(|p| p["format"] == "tsv"));
    for taxon in [9606, 10090] {
        let response = client
            .post(format!("{url}/enrich-proteins"))
            .json(&json!({"proteins":["ATR","Q13535","ALIAS"],"organism_taxon":taxon}))
            .send()
            .await
            .unwrap();
        assert_eq!(response.status(), StatusCode::OK);
        let result: Value = response.json().await.unwrap();
        assert_eq!(result["input_count"], 1);
        assert_eq!(result["organism_taxon"], taxon);
    }
    assert_eq!(
        calls
            .lock()
            .unwrap()
            .iter()
            .filter(|p| p["format"] == "tsv")
            .count(),
        2
    );
    let response = client
        .post(format!("{url}/enrich-proteins"))
        .json(&json!({"proteins":["TP53","Q00002"],"organism_taxon":9606}))
        .send()
        .await
        .unwrap();
    assert_eq!(response.status(), StatusCode::BAD_REQUEST);
    let error: Value = response.json().await.unwrap();
    assert_eq!(error["error"]["code"], "outside_annotation_universe");
    assert_eq!(error["error"]["identifiers"], json!(["Q00002"]));
    // A cached batch from another species must be validated, never reinterpreted by symbol.
    let mut record = ProteinResult::missing("ATR", None);
    record.found = true;
    record.uniprot_id = Some("MOUSEONLY".into());
    record.organism = Some("Mus musculus".into());
    let id = app
        .cache
        .put_analysis(None, json!({"items":[record]}))
        .await
        .unwrap();
    let response = client
        .post(format!("{url}/enrich-proteins"))
        .json(&json!({"result_id":id,"organism_taxon":9606}))
        .send()
        .await
        .unwrap();
    assert_eq!(response.status(), StatusCode::BAD_REQUEST);
    task.abort();
    upstream_task.abort();
}

#[tokio::test]
#[ignore = "live UniProt integration; downloads complete taxon proteome snapshot"]
async fn live_uniprot_exact_ten_protein_set() {
    let dir = tempfile::tempdir().unwrap();
    let app = App {
        cache: Cache::open(dir.path().join("live.sqlite").to_str().unwrap(), 30).unwrap(),
        upstream: Upstream::new(
            "https://rest.uniprot.org".into(),
            "https://www.ebi.ac.uk/QuickGO/services".into(),
        )
        .unwrap(),
        kegg: protein_tools::kegg::Kegg::new("https://rest.kegg.jp".into()).unwrap(),
        permits: Arc::new(Semaphore::new(4)),
    };
    let response=app.enrich(serde_json::from_value(json!({"proteins":LIVE_SYMBOLS,"organism_taxon":9606,"source":"go_bp","background":"proteome"})).unwrap()).await.unwrap();
    assert_eq!(response.input_count, 10);
    let full: Value = app
        .cache
        .get_analysis(&response.result_id)
        .await
        .unwrap()
        .unwrap();
    assert_eq!(full["input_identifiers"], json!(LIVE_SYMBOLS));
    for (i, id) in LIVE_ACCESSIONS.iter().enumerate() {
        assert_eq!(full["identifier_mapping"][i]["canonical"], json!(id));
    }
    println!(
        "live release={} input_count={} N={} proteome_entries={}",
        full["uniprot_release"],
        response.input_count,
        response.background_count,
        full["universe"]["proteome_entry_count"]
    );
}

async fn proteome_fixture(Query(params): Query<HashMap<String, String>>) -> axum::Json<Value> {
    let taxon: u64 = params["query"]
        .strip_prefix("taxonomy_id:")
        .unwrap()
        .parse()
        .unwrap();
    axum::Json(
        json!({"results":[{"id":"UP000005640","proteomeType":"Reference proteome","taxonomy":{"taxonId":taxon}}]}),
    )
}

#[tokio::test]
async fn discovery_never_guesses_a_proteome_or_uses_descendants() {
    for mode in [
        "missing",
        "ambiguous",
        "descendant",
        "incomplete",
        "failure",
        "malformed",
    ] {
        let (app, _calls, _dir, upstream_task) = setup("").await;
        let (discovery_url, discovery_task) = serve(Router::new().route("/proteomes/search", get(move || async move {
            let mut headers=HeaderMap::new();
            let record=json!({"id":"UP000005640","proteomeType":"Reference proteome","taxonomy":{"taxonId":9606}});
            let records=match mode {
                "missing"=>vec![],
                "ambiguous"=>vec![record,json!({"id":"UP000000001","proteomeType":"Reference proteome","taxonomy":{"taxonId":9606}})],
                "descendant"=>vec![json!({"id":"UP000005640","proteomeType":"Reference proteome","taxonomy":{"taxonId":9999}})],
                _=>vec![record],
            };
            if mode=="incomplete" { headers.insert("link","<https://rest.uniprot.org/proteomes/search?cursor=next>; rel=\"next\"".parse().unwrap()); }
            (if mode=="failure" {StatusCode::SERVICE_UNAVAILABLE} else {StatusCode::OK},headers,
                if mode=="malformed" {"{}".into()} else {json!({"results":records}).to_string()})
        }))).await;
        let mut app = app;
        app.upstream = Upstream::new(discovery_url.clone(), discovery_url).unwrap();
        let error = app
            .enrich(
                serde_json::from_value(json!({"proteins":["P00000"],"organism_taxon":9606}))
                    .unwrap(),
            )
            .await
            .err()
            .unwrap();
        assert_eq!(error.0, StatusCode::BAD_GATEWAY, "{mode}");
        assert!(
            app.cache
                .get_analysis::<Value>("go_bp_proteome_taxon_9606_v2")
                .await
                .unwrap()
                .is_none()
        );
        discovery_task.abort();
        upstream_task.abort();
    }
}
