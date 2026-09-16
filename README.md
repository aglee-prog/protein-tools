# protein-tools

Small Rust OpenAPI tool server for authoritative UniProt protein facts and QuickGO annotations. It performs exact lookup and normalization only, with no biological interpretation or classification.

## Run

```sh
docker compose up --build -d
curl http://127.0.0.1:8091/health
curl http://127.0.0.1:8091/openapi.json
curl http://127.0.0.1:8091/protein-info \
  -H 'Content-Type: application/json' \
  -d '{"proteins":["TP53","EGFR","P31749"],"organism":"Homo sapiens"}'
```

Compose binds only to localhost and persists SQLite in `./data`. No Open WebUI configuration is performed. The generated OpenAPI document uses the operation ID `lookupProteinInfo`.

Send a large protein set in **one API call**, for example the included 70-human-protein request:

```sh
curl --fail-with-body http://127.0.0.1:8091/protein-info \
  -H 'Content-Type: application/json' \
  --data-binary @examples/large-request.json
```

Send the complete list; do not manually split it into multiple tool calls. The service handles bounded processing and caching internally. Repeat the same call to reuse successful cached results.

For development without Docker:

```sh
CACHE_DB=./data/cache.sqlite BIND_ADDR=127.0.0.1:8080 cargo run
```

Configuration:

| Variable | Default | Purpose |
| --- | --- | --- |
| `CACHE_DB` | `/data/cache.sqlite` | SQLite database; parent directories created on startup |
| `CACHE_TTL_DAYS` | `30` | Nonnegative whole days; zero disables cache hits |
| `BIND_ADDR` | `0.0.0.0:8080` | Server listen address |
| `RUST_LOG` | `protein_tools=info` | Structured JSON log filtering |

## API and correctness

`POST /protein-info` takes 1–500 gene symbols or UniProt accessions and an optional scientific organism name. The constant `MAX_REQUEST_PROTEINS = 500` is a public safety limit; empty lists and lists exceeding it return HTTP 400. It returns an array in input order, including failures and duplicates. `GET /health` checks only this service. `GET /openapi.json` is generated from Rust types and operation definitions.

Resolution tries the exact accession first, then `gene_exact`. Invalid-accession HTTP 400/404 responses allow the gene strategy to proceed. Returned identifiers and organism scientific names are checked for exact case-insensitive equality. Exact gene names or explicitly listed gene synonyms may match. A unique reviewed Swiss-Prot entry is preferred; multiple remaining records are ambiguous. Generic search is never used. Candidate sets exceeding one 500-record page fail conservatively rather than selecting from incomplete results.

`found` indicates successful UniProt resolution. When QuickGO fails, protein facts remain available with an explicit `error`, an empty annotation list, and no cache write. Consumers must check `error` before treating annotations as complete. Not-found has `found: false` and `error: null`; failures and ambiguity have an error message.

QuickGO pages are fetched sequentially, up to 100 pages of 200 records. Exceeding the limit fails explicitly without caching a partial result. Only annotations for the exact resolved accession are retained. Duplicate normalized tuples are removed. GO qualifiers are retained alongside GO ID, aspect, and evidence, including negation. References and other raw annotation metadata are not returned. No interpretation or term enrichment is performed.

The final successful normalized result is cached under trimmed, case-normalized query and organism keys. TTL uses write time; reads do not refresh it. On hits the original request query is restored and `cached` is true, without upstream requests. Not-found and incomplete/error results are not cached. SQLite uses WAL, a busy timeout, and a mutex on a shared connection; database operations run on Tokio's blocking pool. Cache errors are logged and lookups remain available. Expired keys are replaced on a subsequent successful lookup; there is no background cleanup worker.

One shared HTTPS client has a 5-second connect timeout and 30-second request timeout. There are no automatic retries. `MAX_CONCURRENT_LOOKUPS = 4` controls the shared permits across all requests and each request processes at most four unique entries at once. Each lookup retains the existing per-protein UniProt requests and bounded QuickGO pagination. Normalized duplicates within a request share one lookup, including unresolved/error results, and retain each original query in the output. Their `cached` flags reflect whether that shared result came from SQLite. Simultaneous cache misses across separate requests may still perform duplicate bounded lookups. Large cold requests can take time; clients should allow for upstream latency.

Each completed request logs a concise summary: `requested_count`, `unique_count`, `cache_hits`, `cache_misses`, `resolved_count`, and `unresolved_count`. Cache counts refer to unique identifiers (misses include invalid identifiers that cannot use the cache); resolution counts refer to output entries. `found` determines resolution, including results with incomplete annotations.

## Checks

```sh
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test
```

Tests use temporary SQLite databases and localhost mock upstream servers, with no live internet dependency. They cover normalization, persistence and expiration, response parsing, exact resolution ordering, invalid-accession fallback, reviewed preference, ambiguity, unrelated-record rejection, batch isolation, GO pagination/deduplication, and cache hits with zero upstream calls. They also exercise public size limits, 70/230/500-entry requests, mixed cached/uncached results, normalized duplicates, out-of-order completion, shared concurrency bounds, and the OpenAPI limit.

Upstream API references: [UniProt query fields](https://www.uniprot.org/help/query-fields) and [QuickGO API](https://www.ebi.ac.uk/QuickGO/api/).
