use axum::{
    Router,
    extract::State,
    http::{HeaderMap, StatusCode, Uri},
    routing::get,
};
use protein_tools::{
    api::{App, router},
    cache::Cache,
    kegg::{Kegg, KeggRequest, KeggStatus},
    upstream::Upstream,
};
use serde_json::{Value, json};
use std::{
    collections::{BTreeMap, VecDeque},
    sync::{Arc, Mutex},
    time::{Duration, Instant},
};
use tokio::sync::Semaphore;

#[derive(Clone)]
struct Reply(u16, &'static str, &'static str);
#[derive(Default)]
struct Mock {
    replies: BTreeMap<String, VecDeque<Reply>>,
    calls: Vec<(String, Instant)>,
    taxonomy_records: Option<BTreeMap<String, Value>>,
    queries: Vec<String>,
}
async fn respond(
    State(mock): State<Arc<Mutex<Mock>>>,
    uri: Uri,
) -> (StatusCode, HeaderMap, String) {
    let mut mock = mock.lock().unwrap();
    let path = uri.path().to_owned();
    mock.calls.push((path.clone(), Instant::now()));
    if path == "/uniprotkb/search" && mock.taxonomy_records.is_some() {
        let url = reqwest::Url::parse(&format!("http://localhost{uri}")).unwrap();
        let params: BTreeMap<_, _> = url.query_pairs().into_owned().collect();
        assert_eq!(params["fields"], "accession,organism_id,organism_name");
        assert_eq!(params["format"], "json");
        let query = params["query"].clone();
        let results: Vec<_> = query
            .split('"')
            .enumerate()
            .filter(|(i, _)| i % 2 == 1)
            .filter_map(|(_, id)| mock.taxonomy_records.as_ref().unwrap().get(id).cloned())
            .collect();
        mock.queries.push(query);
        return (
            StatusCode::OK,
            HeaderMap::new(),
            json!({"results": results}).to_string(),
        );
    }

    let reply = mock
        .replies
        .get_mut(&path)
        .and_then(VecDeque::pop_front)
        .unwrap_or(Reply(400, "unexpected request", ""));
    let mut headers = HeaderMap::new();
    if !reply.2.is_empty() {
        headers.insert("retry-after", reply.2.parse().unwrap());
    }
    (
        StatusCode::from_u16(reply.0).unwrap(),
        headers,
        reply.1.to_owned(),
    )
}
struct Fixture {
    url: String,
    mock: Arc<Mutex<Mock>>,
    task: tokio::task::JoinHandle<()>,
    dir: tempfile::TempDir,
    cache: Cache,
    provider: Kegg,
    upstream: Upstream,
}
impl Fixture {
    async fn new() -> Self {
        let mock = Arc::new(Mutex::new(Mock::default()));
        let listener = tokio::net::TcpListener::bind("127.0.0.1:0").await.unwrap();
        let url = format!("http://{}", listener.local_addr().unwrap());
        let app = Router::new()
            .route("/{*path}", get(respond))
            .with_state(mock.clone());
        let task = tokio::spawn(async move {
            axum::serve(listener, app).await.unwrap();
        });
        let dir = tempfile::tempdir().unwrap();
        let cache = Cache::open(dir.path().join("cache.sqlite").to_str().unwrap(), 30).unwrap();
        // These tests exercise the existing pathway stages with organism resolution cached.
        for id in (0..230)
            .map(|n| format!("P{n:05}"))
            .chain(["P04637", "P00533", "P12345", "Q99999"].map(str::to_owned))
        {
            cache
                .put_kegg_value(
                    "",
                    "uniprot-taxonomy-v1",
                    &id,
                    Some(protein_tools::uniprot::ProteinTaxonomy {
                        ncbi_taxonomy_id: Some(9606),
                        organism_name: Some("Homo sapiens".into()),
                    }),
                )
                .await
                .unwrap();
        }
        cache
            .put_kegg("", "taxonomy-codes-v1", "9606", vec!["hsa".into()])
            .await
            .unwrap();
        Self {
            upstream: Upstream::new(url.clone(), url.clone()).unwrap(),
            provider: Kegg::new(url.clone()).unwrap(),
            url,
            mock,
            task,
            dir,
            cache,
        }
    }
    async fn cold() -> Self {
        let f = Self::new().await;
        rusqlite::Connection::open(f.dir.path().join("cache.sqlite"))
            .unwrap()
            .execute("DELETE FROM kegg_cache", [])
            .unwrap();
        f.mock.lock().unwrap().taxonomy_records = Some(BTreeMap::new());
        f
    }
    fn taxonomy(&self, accession: &str, taxid: Option<u64>, name: Option<&str>) {
        self.mock.lock().unwrap().taxonomy_records.as_mut().unwrap().insert(accession.into(), json!({
            "primaryAccession": accession, "organism": {"taxonId": taxid, "scientificName": name}
        }));
    }
    fn reply(&self, path: &str, replies: Vec<Reply>) {
        self.mock
            .lock()
            .unwrap()
            .replies
            .insert(path.into(), replies.into());
    }
    fn calls(&self) -> Vec<String> {
        self.mock
            .lock()
            .unwrap()
            .calls
            .iter()
            .map(|(p, _)| p.clone())
            .collect()
    }
}
impl Drop for Fixture {
    fn drop(&mut self) {
        self.task.abort();
    }
}
fn request(proteins: &[&str]) -> KeggRequest {
    KeggRequest {
        proteins: proteins.iter().map(|s| s.to_string()).collect(),
        pathway_id: None,
        min_proteins: None,
    }
}

#[tokio::test]
async fn batching_persistence_empty_results_and_aggregation() {
    let f = Fixture::new().await;
    f.reply(
        "/conv/hsa/up:P00533+up:P04637+up:P12345+up:Q99999",
        vec![Reply(
            200,
            "up:P00533\thsa:1956\nup:P04637\thsa:7157\nup:P04637\thsa:1\nup:P12345\thsa:2\n",
            "",
        )],
    );
    f.reply("/link/pathway/hsa:1+hsa:1956+hsa:2+hsa:7157", vec![Reply(200, "hsa:1\tpath:hsa04115\nhsa:7157\tpath:hsa04115\nhsa:1956\tpath:hsa04115\nhsa:1956\tpath:hsa05200\n", "")]);
    f.reply(
        "/list/path:hsa04115+path:hsa05200",
        vec![Reply(
            200,
            "hsa04115\tp53 signaling pathway\npath:hsa05200\tPathways in cancer\n",
            "",
        )],
    );
    let inputs = [" p04637 ", "P00533", "P04637", "Q99999", "P12345"];
    let result = f
        .provider
        .lookup(&f.cache, &f.upstream, request(&inputs))
        .await;
    assert!(result.complete);
    assert_eq!(result.proteins.len(), 5);
    assert_eq!(result.proteins[0].query, " p04637 ");
    assert_eq!(result.proteins[0].kegg_gene_ids.len(), 2);
    assert!(result.proteins.iter().all(|p| !p.cached));
    assert!(matches!(result.proteins[3].status, KeggStatus::NoMapping));
    assert!(matches!(result.proteins[4].status, KeggStatus::NoPathways));
    assert_eq!(result.pathways[0].protein_count, 2);
    assert_eq!(result.pathways[0].uniprot_ids, ["P00533", "P04637"]);
    assert_eq!(f.calls().len(), 3);
    // Reopen the database and construct a fresh provider with an unavailable upstream.
    f.task.abort();
    let cache = Cache::open(f.dir.path().join("cache.sqlite").to_str().unwrap(), 30).unwrap();
    let provider = Kegg::new(f.url.clone()).unwrap();
    let mut req = request(&inputs);
    req.min_proteins = Some(2);
    let shared = provider.lookup(&cache, &f.upstream, req).await;
    assert!(shared.complete);
    assert!(shared.proteins.iter().all(|p| p.cached));
    assert_eq!(shared.pathways.len(), 1);
    let mut req = request(&inputs);
    req.pathway_id = Some("hsa05200".into());
    let membership = provider.lookup(&cache, &f.upstream, req).await;
    assert_eq!(membership.pathways[0].uniprot_ids, ["P00533"]);
    assert_eq!(f.calls().len(), 3);
}

#[tokio::test]
async fn partial_names_retry_and_no_cache_poisoning() {
    let f = Fixture::new().await;
    f.reply(
        "/conv/hsa/up:P04637",
        vec![
            Reply(429, "busy", "1"),
            Reply(200, "up:P04637\thsa:7157", ""),
        ],
    );
    f.reply(
        "/link/pathway/hsa:7157",
        vec![
            Reply(503, "unavailable", "0"),
            Reply(200, "hsa:7157\tpath:hsa04115\nhsa:7157\tpath:hsa05200", ""),
        ],
    );
    f.reply(
        "/list/path:hsa04115+path:hsa05200",
        vec![Reply(200, "hsa04115\tp53", "")],
    );
    let result = f
        .provider
        .lookup(&f.cache, &f.upstream, request(&["P04637", "TP53"]))
        .await;
    assert!(!result.complete);
    assert_eq!(result.proteins[0].kegg_gene_ids, ["hsa:7157"]);
    assert_eq!(result.proteins[0].pathways.len(), 2);
    assert!(result.proteins[0].pathways[1].pathway_name.is_none());
    assert!(result.proteins.iter().all(|p| p.error.is_some()));
    assert!(
        f.cache
            .get_kegg("hsa", "names-v1", "hsa05200")
            .await
            .unwrap()
            .is_none()
    );
    f.reply(
        "/list/path:hsa05200",
        vec![Reply(200, "hsa05200\tCancer", "")],
    );
    let recovered = f
        .provider
        .lookup(&f.cache, &f.upstream, request(&["P04637"]))
        .await;
    assert!(recovered.complete);
    assert_eq!(f.calls().len(), 6);
    assert_eq!(f.calls().last().unwrap(), "/list/path:hsa05200");
    let calls = &f.mock.lock().unwrap().calls;
    assert!(calls[1].1.duration_since(calls[0].1) >= Duration::from_millis(990));
}

#[tokio::test]
async fn malformed_and_failed_responses_are_not_cached_and_expired_data_survives() {
    let f = Fixture::new().await;
    f.cache
        .put_kegg("hsa", "conversion-v1", "P04637", vec!["hsa:7157".into()])
        .await
        .unwrap();
    let db = rusqlite::Connection::open(f.dir.path().join("cache.sqlite")).unwrap();
    db.execute(
        "UPDATE kegg_cache SET saved=1 WHERE stage='conversion-v1'",
        [],
    )
    .unwrap();
    f.reply(
        "/conv/hsa/up:P04637",
        vec![
            Reply(200, "<html>outage</html>", ""),
            Reply(503, "down", "0"),
            Reply(503, "down", "0"),
            Reply(503, "down", "0"),
        ],
    );
    for _ in 0..2 {
        let result = f
            .provider
            .lookup(&f.cache, &f.upstream, request(&["P04637"]))
            .await;
        assert!(!result.complete);
        assert!(result.proteins[0].error.is_some());
        assert!(
            f.cache
                .get_kegg("hsa", "conversion-v1", "P04637")
                .await
                .unwrap()
                .is_none()
        );
    }
    let stored: String = db
        .query_row(
            "SELECT result FROM kegg_cache WHERE stage='conversion-v1'",
            [],
            |row| row.get(0),
        )
        .unwrap();
    assert_eq!(stored, "[\"hsa:7157\"]");
    assert_eq!(f.calls().len(), 4);
}

#[tokio::test]
async fn large_batch_chunk_failure_is_isolated() {
    let f = Fixture::new().await;
    let ids: Vec<String> = (0..230).map(|n| format!("P{n:05}")).collect();
    for (i, chunk) in ids.chunks(100).enumerate() {
        let path = format!(
            "/conv/hsa/{}",
            chunk
                .iter()
                .map(|id| format!("up:{id}"))
                .collect::<Vec<_>>()
                .join("+")
        );
        f.reply(
            &path,
            vec![if i == 1 {
                Reply(400, "bad batch", "")
            } else {
                Reply(200, "", "")
            }],
        );
    }
    let result = f
        .provider
        .lookup(
            &f.cache,
            &f.upstream,
            KeggRequest {
                proteins: ids,
                min_proteins: None,
                pathway_id: None,
            },
        )
        .await;
    assert_eq!(f.calls().len(), 3);
    assert_eq!(
        result.proteins.iter().filter(|p| p.error.is_some()).count(),
        100
    );
    assert_eq!(
        result
            .proteins
            .iter()
            .filter(|p| matches!(p.status, KeggStatus::NoMapping))
            .count(),
        130
    );
    assert!(!result.complete);
}

#[tokio::test]
async fn rate_limit_is_global_across_separate_providers_and_concurrent_calls() {
    let f = Fixture::new().await;
    f.reply("/conv/hsa/up:P04637", vec![Reply(200, "", "")]);
    f.reply("/conv/hsa/up:P00533", vec![Reply(200, "", "")]);
    let other = Kegg::new(f.url.clone()).unwrap();
    let (a, b) = tokio::join!(
        f.provider
            .lookup(&f.cache, &f.upstream, request(&["P04637"])),
        other.lookup(&f.cache, &f.upstream, request(&["P00533"]))
    );
    assert!(a.complete && b.complete);
    let calls = &f.mock.lock().unwrap().calls;
    assert_eq!(calls.len(), 2);
    // Allow small localhost scheduling variance between send time and handler arrival.
    assert!(calls[1].1.duration_since(calls[0].1) >= Duration::from_millis(480));
}

#[tokio::test]
async fn public_operation_schema_validation_and_membership() {
    use utoipa::OpenApi;
    let schema = protein_tools::api::ApiDoc::openapi().to_json().unwrap();
    let schema: Value = serde_json::from_str(&schema).unwrap();
    assert_eq!(
        schema["paths"]["/kegg-pathways"]["post"]["operationId"],
        "lookupKeggPathways"
    );
    assert_eq!(
        schema["components"]["schemas"]["KeggRequest"]["properties"]["proteins"]["maxItems"],
        500
    );
    let f = Fixture::new().await;
    f.reply("/conv/hsa/up:P04637", vec![Reply(200, "", "")]);
    let app = App {
        cache: f.cache.clone(),
        kegg: f.provider.clone(),
        upstream: Upstream::new(f.url.clone(), f.url.clone()).unwrap(),
        permits: Arc::new(Semaphore::new(4)),
    };
    let listener = tokio::net::TcpListener::bind("127.0.0.1:0").await.unwrap();
    let url = format!("http://{}/kegg-pathways", listener.local_addr().unwrap());
    let task = tokio::spawn(async move {
        axum::serve(listener, router(app)).await.unwrap();
    });
    let client = reqwest::Client::new();
    for input in [
        json!({"proteins":[]}),
        json!({"proteins":vec!["P04637";501]}),
        json!({"proteins":["P04637"],"min_proteins":0}),
        json!({"proteins":["P04637"],"pathway_id":"mmu:04115"}),
    ] {
        assert_eq!(
            client
                .post(&url)
                .json(&input)
                .send()
                .await
                .unwrap()
                .status(),
            400
        );
    }
    assert!(f.calls().is_empty());
    let response: Value = client
        .post(&url)
        .json(&json!({"proteins":["P04637"],"pathway_id":"path:HSA04115"}))
        .send()
        .await
        .unwrap()
        .json()
        .await
        .unwrap();
    assert_eq!(response["complete"], true);
    assert_eq!(response["proteins"][0]["status"], "no_mapping");
    task.abort();
}

#[tokio::test]
async fn pathway_names_use_ten_entry_batches_and_reuse_cached_links() {
    let f = Fixture::new().await;
    f.cache
        .put_kegg("hsa", "conversion-v1", "P04637", vec!["hsa:7157".into()])
        .await
        .unwrap();
    let paths: Vec<String> = (1..=11).map(|n| format!("hsa{n:05}")).collect();
    f.cache
        .put_kegg("hsa", "links-v1", "hsa:7157", paths.clone())
        .await
        .unwrap();
    let first = format!(
        "/list/{}",
        paths[..10]
            .iter()
            .map(|id| format!("path:{id}"))
            .collect::<Vec<_>>()
            .join("+")
    );
    f.reply(&first, vec![Reply(200, "hsa00001\tOne\nhsa00002\tTwo\nhsa00003\tThree\nhsa00004\tFour\nhsa00005\tFive\nhsa00006\tSix\nhsa00007\tSeven\nhsa00008\tEight\nhsa00009\tNine\nhsa00010\tTen", "")]);
    f.reply(
        "/list/path:hsa00011",
        vec![Reply(200, "hsa00011\tEleven", "")],
    );
    let result = f
        .provider
        .lookup(&f.cache, &f.upstream, request(&["P04637"]))
        .await;
    assert!(result.complete);
    assert_eq!(result.pathways.len(), 11);
    assert_eq!(f.calls().len(), 2);
    assert!(f.calls().iter().all(|p| p.starts_with("/list/")));
    assert!(!result.proteins[0].cached);
    let result = f
        .provider
        .lookup(&f.cache, &f.upstream, request(&["P04637"]))
        .await;
    assert!(result.proteins[0].cached);
    assert_eq!(f.calls().len(), 2);
}

#[tokio::test]
async fn unavailable_network_keeps_successful_cached_proteins() {
    let f = Fixture::new().await;
    f.cache
        .put_kegg("hsa", "conversion-v1", "P04637", vec![])
        .await
        .unwrap();
    f.task.abort();
    // Wait for the mock listener to be dropped, making connection refusal deterministic.
    tokio::task::yield_now().await;
    let result = f
        .provider
        .lookup(&f.cache, &f.upstream, request(&["P04637", "P00533"]))
        .await;
    assert!(!result.complete);
    assert!(result.proteins[0].cached);
    assert!(matches!(result.proteins[0].status, KeggStatus::NoMapping));
    assert!(
        result.proteins[1]
            .error
            .as_ref()
            .unwrap()
            .contains("network")
    );
    assert!(
        f.cache
            .get_kegg("hsa", "conversion-v1", "P00533")
            .await
            .unwrap()
            .is_none()
    );
}

#[tokio::test]
async fn mixed_organisms_resolve_dynamically_and_persist_across_restart() {
    let f = Fixture::cold().await;
    for (id, taxid, name) in [
        ("P04637", 9606, "Homo sapiens"),
        ("P00533", 9606, "Homo sapiens"),
        ("P02340", 10090, "Mus musculus"),
        ("P10361", 10116, "Rattus norvegicus"),
        ("P12345", 999999, "Unsupported organism"),
        ("P54321", 562, "Escherichia coli"),
    ] {
        f.taxonomy(id, Some(taxid), Some(name));
    }
    f.taxonomy("Q99998", None, Some("Unknown taxonomy"));
    f.reply(
        "/link/genome/taxid:10090+taxid:10116+taxid:562+taxid:9606+taxid:999999",
        vec![Reply(
            200,
            "taxid:9606\tgn:hsa\ntaxid:10090\tgn:mmu\ntaxid:10116\tgn:rno\ntaxid:562\tgn:eco",
            "",
        )],
    );
    f.reply(
        "/conv/hsa/up:P00533+up:P04637",
        vec![Reply(200, "up:P04637\thsa:7157\nup:P00533\thsa:1956", "")],
    );
    f.reply(
        "/link/pathway/hsa:1956+hsa:7157",
        vec![Reply(
            200,
            "hsa:7157\tpath:hsa04115\nhsa:1956\tpath:hsa04115",
            "",
        )],
    );
    f.reply(
        "/list/path:hsa04115",
        vec![Reply(
            200,
            "hsa04115\tp53 signaling pathway - Homo sapiens",
            "",
        )],
    );
    f.reply(
        "/conv/mmu/up:P02340",
        vec![Reply(200, "up:P02340\tmmu:22059", "")],
    );
    f.reply(
        "/link/pathway/mmu:22059",
        vec![Reply(200, "mmu:22059\tpath:mmu04115", "")],
    );
    f.reply(
        "/list/path:mmu04115",
        vec![Reply(
            200,
            "mmu04115\tp53 signaling pathway - Mus musculus",
            "",
        )],
    );
    f.reply(
        "/conv/rno/up:P10361",
        vec![Reply(200, "up:P10361\trno:24842", "")],
    );
    f.reply(
        "/link/pathway/rno:24842",
        vec![Reply(200, "rno:24842\tpath:rno04115", "")],
    );
    f.reply(
        "/list/path:rno04115",
        vec![Reply(
            200,
            "rno04115\tp53 signaling pathway - Rattus norvegicus",
            "",
        )],
    );
    // A fourth organism and a non-numeric locus tag demonstrate generality.
    f.reply(
        "/conv/eco/up:P54321",
        vec![Reply(200, "up:P54321\teco:b0001", "")],
    );
    f.reply(
        "/link/pathway/eco:b0001",
        vec![Reply(200, "eco:b0001\tpath:eco00010", "")],
    );
    f.reply(
        "/list/path:eco00010",
        vec![Reply(200, "eco00010\tGlycolysis", "")],
    );
    let inputs = [
        "P04637", "P02340", "P10361", "P12345", "Q99999", "Q99998", "P00533", "P54321", " p04637 ",
    ];
    let result = f
        .provider
        .lookup(&f.cache, &f.upstream, request(&inputs))
        .await;
    assert!(result.complete, "{result:?}");
    for (i, taxid, code, name) in [
        (0, 9606, "hsa", "Homo sapiens"),
        (1, 10090, "mmu", "Mus musculus"),
        (2, 10116, "rno", "Rattus norvegicus"),
    ] {
        let protein = &result.proteins[i];
        assert_eq!(protein.ncbi_taxonomy_id, Some(taxid));
        assert_eq!(protein.organism_name.as_deref(), Some(name));
        assert_eq!(protein.kegg_organism_codes, [code]);
        assert!(matches!(protein.status, KeggStatus::Mapped));
        assert!(protein.pathways[0].pathway_id.starts_with(code));
    }
    assert!(matches!(
        result.proteins[3].status,
        KeggStatus::UnsupportedOrganism
    ));
    assert!(matches!(result.proteins[4].status, KeggStatus::NoMapping));
    assert!(matches!(
        result.proteins[5].status,
        KeggStatus::UnsupportedOrganism
    ));
    assert_eq!(result.proteins[3].ncbi_taxonomy_id, Some(999999));
    assert_eq!(result.proteins[7].kegg_gene_ids, ["eco:b0001"]);
    assert_eq!(result.pathways.len(), 4); // Same numbered pathways remain organism-specific.
    assert_eq!(
        result
            .pathways
            .iter()
            .find(|p| p.pathway_id == "hsa04115")
            .unwrap()
            .protein_count,
        2
    );
    assert_eq!(
        f.calls()
            .iter()
            .filter(|p| *p == "/uniprotkb/search")
            .count(),
        1
    );
    assert_eq!(f.calls().len(), 14);
    assert_eq!(
        f.mock.lock().unwrap().queries[0]
            .matches("accession:")
            .count(),
        8
    );
    assert!(result.proteins.iter().all(|p| !p.cached));
    f.task.abort();
    let cache = Cache::open(f.dir.path().join("cache.sqlite").to_str().unwrap(), 30).unwrap();
    let provider = Kegg::new(f.url.clone()).unwrap();
    // Separate human, mouse, rat, and unsupported/unknown requests all work offline.
    for id in inputs {
        let result = provider.lookup(&cache, &f.upstream, request(&[id])).await;
        assert!(result.complete);
        assert!(result.proteins[0].cached);
    }
    let mut filtered = request(&["P04637", "P02340", "P10361"]);
    filtered.pathway_id = Some("mmu04115".into());
    let result = provider.lookup(&cache, &f.upstream, filtered).await;
    assert_eq!(result.pathways.len(), 1);
    assert_eq!(result.pathways[0].uniprot_ids, ["P02340"]);
    assert_eq!(f.calls().len(), 14);
}

#[tokio::test]
async fn new_accessions_reuse_taxonomy_code_cache_and_resolution_is_batched() {
    let f = Fixture::cold().await;
    let ids: Vec<_> = (0..101).map(|n| format!("P{n:05}")).collect();
    for id in &ids {
        f.taxonomy(id, Some(123456), Some("Not in KEGG"));
    }
    f.reply("/link/genome/taxid:123456", vec![Reply(200, "", "")]);
    let result = f
        .provider
        .lookup(
            &f.cache,
            &f.upstream,
            KeggRequest {
                proteins: ids,
                pathway_id: None,
                min_proteins: None,
            },
        )
        .await;
    assert!(result.complete);
    assert!(
        result
            .proteins
            .iter()
            .all(|p| matches!(p.status, KeggStatus::UnsupportedOrganism))
    );
    let queries = f.mock.lock().unwrap().queries.clone();
    assert_eq!(queries.len(), 2);
    assert_eq!(queries[0].matches("accession:").count(), 100);
    assert_eq!(queries[1].matches("accession:").count(), 1);
    f.taxonomy("Q99999", Some(123456), Some("Not in KEGG"));
    let result = f
        .provider
        .lookup(&f.cache, &f.upstream, request(&["Q99999"]))
        .await;
    assert!(result.complete);
    assert!(matches!(
        result.proteins[0].status,
        KeggStatus::UnsupportedOrganism
    ));
    assert_eq!(
        f.calls()
            .iter()
            .filter(|p| p.starts_with("/link/genome/"))
            .count(),
        1
    );
    assert_eq!(
        f.calls()
            .iter()
            .filter(|p| *p == "/uniprotkb/search")
            .count(),
        3
    );
}

#[tokio::test]
async fn organism_resolution_failures_are_not_cached_or_called_unsupported() {
    let f = Fixture::cold().await;
    f.mock.lock().unwrap().taxonomy_records = None;
    f.reply("/uniprotkb/search", vec![Reply(503, "unavailable", ""), Reply(200, "{broken", ""), Reply(200, "{\"results\":[{\"primaryAccession\":\"P02340\",\"organism\":{\"taxonId\":10090,\"scientificName\":\"Mus musculus\"}}]}", "")]);
    for _ in 0..2 {
        let result = f
            .provider
            .lookup(&f.cache, &f.upstream, request(&["P02340"]))
            .await;
        assert!(!result.complete);
        assert!(matches!(result.proteins[0].status, KeggStatus::Error));
        assert!(
            f.cache
                .get_kegg_value::<Option<protein_tools::uniprot::ProteinTaxonomy>>(
                    "",
                    "uniprot-taxonomy-v1",
                    "P02340"
                )
                .await
                .unwrap()
                .is_none()
        );
    }
    f.reply(
        "/link/genome/taxid:10090",
        vec![
            Reply(200, "<html>bad response</html>", ""),
            Reply(200, "taxid:10090\tgn:mmu", ""),
        ],
    );
    let result = f
        .provider
        .lookup(&f.cache, &f.upstream, request(&["P02340"]))
        .await;
    assert!(!result.complete);
    assert_eq!(result.proteins[0].ncbi_taxonomy_id, Some(10090));
    assert!(
        f.cache
            .get_kegg("", "taxonomy-codes-v1", "10090")
            .await
            .unwrap()
            .is_none()
    );
    f.reply("/conv/mmu/up:P02340", vec![Reply(200, "", "")]);
    let result = f
        .provider
        .lookup(&f.cache, &f.upstream, request(&["P02340"]))
        .await;
    assert!(result.complete);
    assert_eq!(result.proteins[0].kegg_organism_codes, ["mmu"]);
    assert_eq!(
        f.calls()
            .iter()
            .filter(|p| *p == "/uniprotkb/search")
            .count(),
        3
    );
}
