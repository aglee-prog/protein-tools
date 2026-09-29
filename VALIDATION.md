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

## Stage 1 GO BP enrichment — 2026-09-28

- `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`, and `git diff --check`: passed.
- Full `cargo test`: passed (6 unit tests, 3 enrichment integration tests, 14 KEGG integration tests, 11 service integration tests). After adding one more cached-species regression and extra pagination failure cases, `cargo test --test enrichment` passed all 4 enrichment integration tests; formatting and Clippy were rerun successfully.
- New deterministic checks cover exact hypergeometric upper tails against combinatorial values, BH adjustment across zero-hit terms, monotonicity, population boundaries, unannotated proteins, root exclusion, duplicate memberships, compact responses, persisted full results, aliases, restart and expiration, OpenAPI registration, and cached protein batches.
- Localhost upstream mocks verify complete two-page loading, comma-containing Link URLs, failed/truncated/malformed pages, release changes, repeated cursors, duplicate background entries, missing population headers and empty annotation catalogs. Invalid snapshots never enter the cache. Cached unresolved or nonhuman protein sets cannot be reinterpreted as human gene symbols.
- A live official UniProt request confirmed the TSV `accession,go_p` format, background query, pagination Link and population/release headers. A complete live proteome download and end-to-end live enrichment were not run. Automated tests use local fixtures and no live upstream APIs; socket-restricted sandbox execution requires the approved test command outside the sandbox.

## Stage 2 selective enrichment inspection — 2026-09-28

- `cargo test`: 39 offline tests passed (6 unit, 8 enrichment, 14 KEGG, 11 service); the opt-in live test was run separately and passed. After adding expiry assertions, all 8 offline enrichment tests passed again.
- `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`, and `git diff --check`: passed.
- Stage 1 fixtures now verify default and capped term limits, stable FDR/p-value/GO-ID ordering, offsets and end pages, single-term statistics, selected hit membership, aliases and original input order, zero-hit terms, invalid/missing/expired/wrong-type IDs, SQLite reopen, unchanged small/protein retrieval, large-result blocking, and OpenAPI registration. Inspection performs no upstream requests.
- Live command: `cargo test --test enrichment live_uniprot_exact_ten_protein_set -- --ignored --nocapture`. This uses a temporary HTTP service running the changed code, a fresh temporary SQLite cache, and official UniProt APIs; it does not deploy or modify an existing service.
- Exact input: TP53, ATM, ATR, BRCA1, BRCA2, CHEK1, CHEK2, RAD51, MDM2, CDKN1A; taxon 9606. Completed in 90.98 seconds with 11,021 tested terms.
- HTTP workflow response sizes: enrichment summary **1,641 bytes**; first 20 terms **4,404 bytes**; selected term proteins **556 bytes**; protected generic retrieval **275 bytes**, with `truncated=true` and no full payload.
- Selected GO:0006974, “DNA damage response”: nine canonical hits, with original mappings for TP53, ATM, ATR, BRCA1, BRCA2, CHEK1, CHEK2, RAD51 and CDKN1A. The live test checks these mappings against the established Stage 1 expected accessions. No complete enrichment dataset crosses the HTTP boundary in this workflow.
- Persistence format and statistical calculations are unchanged. Selective operations deserialize the saved snapshot inside the service; the new bounds protect outgoing responses rather than changing server-side storage.
