# Validation — 2026-09-16

- `cargo fmt --check`: passed.
- `cargo clippy --all-targets --all-features -- -D warnings`: passed.
- `cargo test`: all 5 integration tests passed. Mock HTTP servers require localhost socket permission; no external API access is used by tests.
- Live checks ran the Rust binary on `127.0.0.1:8092` with `CACHE_DB=./data/workstation-check.sqlite`.
- `GET /health`: HTTP 200, `{"status":"ok"}`.
- `GET /openapi.json`: HTTP 200, OpenAPI 3.1.0; generated request schema documents gene symbols/accessions, optional organism, and 1–50 inputs.
- First TP53 + Homo sapiens lookup: `found=true`, `uniprot_id=P04637`, `gene=TP53`, `organism=Homo sapiens`, protein name `Cellular tumor antigen p53`, `cached=false`, `error=null`, 213 normalized GO annotations.
- Identical repeat: `cached=true`, with a cache-hit event and no additional UniProt/QuickGO events for TP53/P04637.
- Batch TP53, EGFR, P31749: all resolved without errors to P04637/TP53, P00533/EGFR, P31749/AKT1 respectively, all Homo sapiens. TP53 was cached. Annotation counts were 213, 138, 201.

Docker image build and compose startup remain unverified: the workstation user lacks permission to access `/var/run/docker.sock`, including outside the execution sandbox. Passwordless sudo is unavailable (`sudo: a password is required`). No system permissions or Docker configuration were changed. Run `docker compose up --build -d` from an account with Docker access to complete container verification.
