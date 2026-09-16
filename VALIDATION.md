# Validation — 2026-09-16

- `cargo fmt --check`: passed.
- `cargo clippy --all-targets --all-features -- -D warnings`: passed.
- `cargo test`: all 8 integration tests passed. Mock HTTP servers require localhost socket permission; tests passed outside the socket-restricted sandbox and use no live upstream APIs.
- Existing resolver/cache/QuickGO tests remain. Added public HTTP size checks (0 and 501 rejected; 70, 200, 230, and 500 accepted), mixed cache batches including 200 cached + 30 fresh, repeat requests with zero upstream calls, normalized duplicate success/failure sharing, OpenAPI limits, forced out-of-order completion, and a shared four-lookup bound across concurrent requests.

## Live verification

Ran the updated Rust binary on `127.0.0.1:8093` with a separate initially empty `CACHE_DB=./data/large-request-check.sqlite`.

- `GET /health`: HTTP 200, `{"status":"ok"}`.
- `GET /openapi.json`: HTTP 200; schema has `minItems: 1`, `maxItems: 500`, and operation wording instructs clients to send the complete list in one request.
- One POST using `examples/large-request.json` (70 known human gene symbols): HTTP 200 in 80.99 seconds, 70 entries in exact input order, 68 resolved, 0 cache hits. MET and CDKN2A returned existing ambiguous-exact-match errors. Their failures did not discard the other 68 results.
- Identical repeat: HTTP 200 in 2.25 seconds, same order, 68 cache hits. Upstream lookup/page events increased only from 282 to 286: accession and gene_exact for each of the two ambiguous identifiers. No cached protein triggered UniProt or QuickGO calls, and no QuickGO pages were requested on the repeat. Unresolved results remain uncached under the existing policy.
- Additional single POST with `NOT_A_REAL_PROTEIN` inserted as the second item: HTTP 200 in 2.70 seconds, 71 entries in exact input order, 68 successful cache hits. The inserted entry had `found: false`, `error: null`; the two ambiguity errors remained isolated.
- Structured request summaries agreed with the returned counts. Response files and logs are retained in ignored `data/large-*.json` and `data/large-request-check.log`.

## Container limitation

`docker compose up --build -d` was attempted both inside and outside the sandbox, but the workstation user cannot access `/var/run/docker.sock`. `sudo -n docker compose up --build -d` also failed: `sudo: a password is required`. Therefore Docker image build and Compose startup remain unverified; live checks above used the native binary. No system permissions or Docker configuration were changed. Run `docker compose up --build -d` from an account with Docker access to complete container verification.

Large cold requests remain sensitive to upstream latency and client timeouts. Deduplication is per request; simultaneous separate requests may duplicate bounded work. No resolver, annotation, or cache policy was changed.
