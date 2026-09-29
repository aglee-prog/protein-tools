# protein-tools

Small Rust OpenAPI tool server for authoritative UniProt protein facts and QuickGO annotations. It performs exact lookup and normalization, plus deterministic organism-scoped GO Biological Process overrepresentation analysis.

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

For broad functional analysis, request compact identity, function, and Biological Process data:

```json
{
  "proteins": ["TP53", "EGFR", "AKT1", "..."],
  "organism": "Homo sapiens",
  "include": ["identity", "function", "go.biological_process"],
  "max_go_terms_per_protein": 10
}
```

For identity and function only:

```json
{
  "proteins": ["TP53", "EGFR"],
  "organism": "Homo sapiens",
  "include": ["identity", "function"]
}
```

For a detailed follow-up:

```json
{
  "proteins": ["TP53"],
  "organism": "Homo sapiens",
  "include": ["identity", "function", "go.biological_process", "go.molecular_function", "go.cellular_component", "go.evidence"],
  "max_go_terms_per_protein": 50
}
```

Detailed requests with evidence are best directed at a small subset of proteins because the response grows with the annotations returned.

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

`include` is optional and accepts only `identity`, `function`, `go.biological_process`, `go.molecular_function`, `go.cellular_component`, and `go.evidence`. The default is identity, function, and Biological Process. Every result always includes `query`, `found`, `cached`, and `error`; other fields appear only when selected and available. `go.evidence` adds available QuickGO evidence codes to the selected GO aspects and does not select an aspect by itself. GO evidence is otherwise omitted. If `max_go_terms_per_protein` is omitted, all matching GO annotations are returned. An explicit value of at least 1 limits the total annotations across selected aspects for each protein. Annotations are sorted by aspect, GO ID, qualifier, and evidence before truncation. This ordering is deterministic and does not rank biological importance. Invalid selection values or a zero limit return HTTP 400.

Resolution tries the exact accession first, then `gene_exact`. Invalid-accession HTTP 400/404 responses allow the gene strategy to proceed. Returned identifiers and organism scientific names are checked for exact case-insensitive equality. Exact gene names or explicitly listed gene synonyms may match. A unique reviewed Swiss-Prot entry is preferred; multiple remaining records are ambiguous. Generic search is never used. Candidate sets exceeding one 500-record page fail conservatively rather than selecting from incomplete results.

`found` indicates successful UniProt resolution. When QuickGO fails, protein facts remain available with an explicit `error`, an empty annotation list, and no cache write. Consumers must check `error` before treating annotations as complete. Not-found has `found: false` and `error: null`; failures and ambiguity have an error message.

QuickGO pages are fetched sequentially, up to 100 pages of 200 records. Exceeding the limit fails explicitly without caching a partial result. Only annotations for the exact resolved accession are retained. Duplicate normalized tuples are removed. GO qualifiers are retained alongside GO ID, aspect, and evidence, including negation. References and other raw annotation metadata are not returned. The protein-info endpoint performs no enrichment; use enrich_proteins for deterministic GO BP analysis.

The final successful normalized result is cached under trimmed, case-normalized query and organism keys. Selection and GO limits are applied after reading the complete record, so they do not create separate cache entries. A fresh identity/function-only request still fetches QuickGO to populate that complete record for later requests. TTL uses write time; reads do not refresh it. On hits the original request query is restored and `cached` is true, without upstream requests. Not-found and incomplete/error results are not cached. SQLite uses WAL, a busy timeout, and a mutex on a shared connection; database operations run on Tokio's blocking pool. Cache errors are logged and lookups remain available. Expired keys are replaced on a subsequent successful lookup; there is no background cleanup worker.

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

## KEGG pathways (mixed organisms)

`POST /kegg-pathways`, OpenAPI operation **`lookupKeggPathways`**, accepts 1–500
canonical UniProt accessions from any organism, including mixed-organism sets.
This repository exposes OpenAPI operations for tool consumers; it does not implement a native MCP transport. Full protein and enrichment results are available through the SQLite `result_id` store described below.
Existing UniProt/QuickGO behavior is unchanged.

```json
{
  "proteins": ["P04637", "P02340", "P10361"],
  "min_proteins": 1
}
```

Omit `min_proteins` to obtain every associated pathway. Use 2 for pathways shared
by at least two distinct supplied accessions, or the number of distinct supplied
accessions for pathways common to the entire set. Add `"pathway_id": "hsa04115"`
to obtain supplied-set membership for that pathway (`path:hsa04115` also works).
Both filters apply only to the `pathways` groups, leaving the per-protein results
available for inspection. Repeating inputs does not increase membership counts.
Pathway IDs retain their organism code: `hsa04115`, `mmu04115`, and `rno04115`
are separate groups, even though their pathway numbers match. Shared-pathway
counts use exact IDs; cross-species orthology or equivalence is not inferred.
These are associations reported by KEGG, not statistical enrichment results.

The structured response contains:

- `complete`, which is false if any input failed. Organism metadata is per protein;
  the former top-level single `organism` field has been removed.
- `proteins`, in original input order: `query`, normalized `uniprot_id`,
  `ncbi_taxonomy_id`, `organism_name`, `kegg_organism_codes`, `kegg_gene_ids`,
  gene/pathway associations with `pathway_id` and `pathway_name`,
  `status`, `cached`, `error`, and `source: "KEGG"`.
- `pathways`, sorted by ID: pathway name, distinct sorted `uniprot_ids`,
  `protein_count`, and source. Only proteins supplied in the request are included.

Organisms are resolved from UniProt metadata using batches of up to 100 exact
accessions (`accession`, `organism_id`, and `organism_name` fields). KEGG's
`link/genome/taxid:<id>+...` endpoint then supplies the organism codes for each
exact NCBI taxonomy ID. There is no hard-coded species table, name matching, or
fallback to a related species/strain. All returned organism codes are preserved
because a taxonomy can link to multiple KEGG genomes. Accessions are grouped by
those codes before conversion and pathway requests. Non-numeric gene locus tags
are supported as well as NCBI GeneIDs.

Statuses are `mapped`, `no_mapping`, `no_pathways`, `unsupported_organism`, and
`error`. An unknown primary accession or a supported organism with no KEGG gene
mapping returns `no_mapping`. A record without a taxonomy ID, or with no KEGG
organism linked to its exact taxonomy, returns `unsupported_organism`. Both are
normal per-protein results (`error: null`) and do not mark the batch incomplete.
Neither proves a lack of biological pathway involvement. Upstream failures
remain errors and are never cached as unsupported organisms. Gene symbols,
secondary accessions, and isoform suffixes are not silently converted to primary
accessions. Trimmed, case-normalized duplicates share work but retain their
original queries in the output. Multiple gene and organism mappings are retained.

Check `complete` and per-protein `error` before interpreting an empty result or
an intersection. Batch failures are isolated to affected identifiers; successful
mappings and links remain available. A failed pathway-name lookup leaves its
known association present with `pathway_name: null` and an explicit error.
KEGG HTTP 4xx other than 429 and malformed responses fail explicitly. KEGG network
errors, 429, and 5xx receive at most three attempts, with backoff and Retry-After
support. UniProt metadata requests use the existing bounded client without retries.
A server cooldown longer than 30 seconds returns an error promptly and blocks
new upstream attempts until the cooldown permits them; cached data remains usable.

KEGG uses the existing `CACHE_DB` SQLite file and `CACHE_TTL_DAYS` (default 30).
The existing `kegg_cache` table stores normalized per-accession UniProt taxonomy
metadata, per-taxonomy KEGG organism codes, and per-identifier conversion, link,
and name results. Organism and stage/version remain part of the keys; taxonomy
resolution uses a separate global namespace. Successful unknown accessions and
unsupported taxonomy results are also cached. Repeated lookups, including after
restart, need no UniProt or KEGG requests while all required entries are fresh.
New accessions with an already cached taxonomy reuse its organism-code mapping.
Existing conversion/link/name cache entries remain compatible. Each stage reads
its cache before batching misses. Successful empty conversions and links are
cached; errors, malformed bodies, and missing names are not. Writes do not
replace other stages, and failures leave older entries untouched. Expired data
is not served as fresh. Cache failures are logged and lookups continue. A
protein's `cached` flag is true only if every stage it needs came from SQLite.
Persistence survives restarts; there is no external cache service or cleanup
worker. Concurrent cold requests can duplicate a miss, but share the rate limit.

The provider uses only the official `https://rest.kegg.jp` API:

- `/link/genome/taxid:<taxonomy>+taxid:<taxonomy>...` (up to 100 distinct taxonomies).
- `/conv/<organism>/up:<accession>+up:<accession>...` (up to 100 inputs per organism).
- `/link/pathway/<organism>:<gene>+...` (up to 100 inputs per organism).
- `/list/pathway/<organism>` (one organism catalog when pathway names are missing;
  only names for requested memberships are retained).

Membership IDs (`path:hsa04010`) and catalog IDs (`hsa04010`) normalize to the
same cache key. Invalid cached names, including empty or null values, are retried
without invalidating taxonomy, conversion, membership, or valid name entries.
Only requested gene mappings and memberships are retrieved. All KEGG provider instances and retries in **one service
process** share a conservative minimum 500 ms request interval (2 requests/sec).
Multiple replicas do not coordinate their limits; run one instance per shared
outbound KEGG traffic budget. Resolution uses exact taxonomy links, so an organism
without such a link is reported as unsupported even if KEGG represents a related
strain or species. The tiny `httpdate` dependency, already present transitively, parses HTTP-date Retry-After
headers. Tests use only localhost mocks and temporary SQLite databases.

See the [official KEGG API manual](https://www.genome.jp/kegg/rest/keggapi.html)
for endpoint formats and the [KEGG API access conditions](https://www.genome.jp/kegg/rest/)
for academic-use terms and the published request limit.

## Compact KEGG comparison

`POST /compare-kegg-pathways`, OpenAPI operation **`compareKeggPathways`**, reuses
`lookupKeggPathways` and its persistent cache. Prefer this operation for comparisons
that need shared pathways without the raw per-protein associations.

Request (only `proteins` is required):

```json
{"proteins":["P31749","Q05030"],"min_proteins":2}
```

`proteins` accepts 2–500 strings; canonical accession validation and normalization
follow the lookup. Invalid entries become per-input errors, allowing valid entries
to contribute. `min_proteins` is an optional integer from 1–500 (null or omitted
means 2). Set it to the distinct accession count for intersection; thresholds above
the input count return no matches. Duplicate inputs do not increase counts.

Exact response shape (all fields are always present):

```text
{
  complete: boolean,
  protein_count: integer,
  matched_pathway_count: integer,
  proteins: [{
    query: string,
    uniprot_id: string,
    status: "mapped" | "no_mapping" | "unsupported_organism" | "no_pathways" | "error",
    cached: boolean,
    error: string | null
  }],
  pathways: [{
    pathway_key: string,
    pathway_name: string | null,
    proteins: [{uniprot_id: string, kegg_gene_id: string, pathway_id: string}],
    protein_count: integer
  }]
}
```

Top-level `protein_count` counts distinct syntactically valid canonical accessions,
including unmapped or failed ones. Statuses retain input order and duplicates.
`complete` follows lookup semantics: false if any input has an error; unmapped or
unsupported accessions are successful negative results. Partial known associations
still contribute. Missing names remain null without discarding membership.

Each pathway's `protein_count` counts distinct accessions. Its `proteins` contains
unique accession/gene/pathway tuples for provenance (multiple genes can represent
one protein). Groups are sorted by key, members lexically by accession/gene/pathway.
Only groups meeting the threshold are returned; status entries have no raw pathways.

Normalization removes the resolved KEGG organism code from each pathway ID and
validates the remaining five digits, preserving leading zeros. Thus `hsa04151`,
`mmu04151`, and `rno04151` share key `04151`; four-character organism codes work too.
Names never establish equivalence. The final ` - organism` catalog suffix is removed
for display; if normalized names differ, the lexically first available name is used.
No cache tables, cache keys, or existing lookup request/response schemas change.


## Stage 1: GO Biological Process enrichment

`POST /enrich-proteins`, OpenAPI operation **`enrich_proteins`**, accepts exactly
one of `proteins` (1–500 gene symbols or UniProt accessions) or `result_id`.
Use `organism_taxon` (an NCBI taxonomy ID) for every new analysis, including
accession inputs. `source="go_bp"` and `background="proteome"` are the supported
defaults. A cached enrichment result supplies its own taxon; a conflicting taxon
is rejected. There is no implicit human default.
```sh
curl --fail-with-body http://127.0.0.1:8091/enrich-proteins \
  -H 'Content-Type: application/json' \
  -d '{"proteins":["TP53","EGFR","AKT1"],"organism_taxon":9606,"source":"go_bp","background":"proteome"}'
```

The compact response contains `result_id`, `organism_taxon`, `method`, `source`, `background`,
`input_count`, `background_count`, `tested_terms`, `significant_terms`,
`fdr_threshold` (0.05), and at most ten `top_terms` with GO ID, name, hit count,
p-value and FDR. Terms are ordered by FDR, p-value, then GO ID. Top terms may
be nonsignificant; use their FDR rather than interpreting their presence as
significance. Protein lists and the full term collection stay outside this response.

Inspect enrichment progressively using the saved `result_id`:

- `GET /results/{result_id}/enrichment/terms?limit=20&offset=0`
  (`get_enrichment_terms`) returns term statistics without membership or evidence.
  Ordering is FDR ascending, then p-value, then GO ID. Default limit is 20;
  larger limits are clamped to 50 and zero is rejected. The response includes
  `total`, effective `limit`, `offset`, `has_more`, and `items`. Terms may be nonsignificant.
- `GET /results/{result_id}/enrichment/terms/{term_id}`
  (`get_enrichment_term`) returns one term's statistics, taxon, method, source,
  background policy and UniProt release.
- `GET /results/{result_id}/enrichment/terms/{term_id}/proteins`
  (`get_enrichment_term_proteins`) returns only that term's canonical input hits
  and original-to-canonical `identifier_mapping`, preserving aliases, repeated
  inputs and input order. `hit_count` counts distinct canonical proteins.
  No extra upstream lookups occur. Use existing protein tools explicitly if evidence is needed.

For example, after enriching the ten DNA damage response genes:
```sh
curl --fail-with-body "http://127.0.0.1:8091/results/$RESULT_ID/enrichment/terms?limit=20"
# Select a GO ID from the returned terms:
curl --fail-with-body "http://127.0.0.1:8091/results/$RESULT_ID/enrichment/terms/$TERM_ID/proteins"
```

`GET /results/{result_id}` (`getCachedResult`) preserves complete small results
and existing protein-batch behavior. Enrichment results exceeding **32 KiB of
serialized JSON** return metadata only: `result_id`, `type="go_bp_enrichment"`,
`item_count`, `truncated=true`, and a message directing clients to the selective
operations. This is explicitly not the full dataset; no partial term list is
passed off as complete. Enrichment inspection should always use the selective tools.
The complete saved result (including universe accessions) stays in SQLite with
its existing TTL and remains selectively inspectable after service restart.
Invalid IDs and wrong result types return HTTP 400; missing/expired results and
absent terms return HTTP 404. No statistics are recalculated during inspection.
`POST /enrich-proteins` with `{"result_id":"result_..."}` still reanalyzes that
protein set using the currently cached background.

`POST /protein-info` retains its existing JSON array response and additionally
returns `X-Result-Id` when storage succeeds. That ID identifies the full unfiltered
protein batch (`proteins` and `items`) and can be supplied to enrichment. Provide `organism_taxon` when using a protein-info result ID. Cached
batches containing unresolved records are rejected; their canonical accessions
are validated against the requested taxon, preserving the original query mapping. If saving a protein batch fails,
the lookup still succeeds without the header; an enrichment cache-write failure
returns HTTP 500 rather than an unusable result ID. Tool clients that do not
expose response headers can supply `proteins` directly.

### Population and annotation policy

The analysis unit is a **distinct primary UniProtKB accession**, representing a
canonical entry, not a gene or an expanded isoform sequence. The background query
is `proteome:<discovered_id> AND organism_id:<organism_taxon>`. The service
discovers the exact taxon’s unique reference proteome using UniProt’s proteomes
API; missing, ambiguous or incomplete discovery fails explicitly. The selected
ID is recorded as `universe.proteome_id`. Entries are deduplicated by accession. Both reviewed and unreviewed
entries are eligible. Multiple entries associated with the same gene remain
separate units; no one-protein-per-gene selection is implied.

Only entries with at least one non-root GO BP association are eligible for this
annotation-conditioned analysis. `N = background_count = universe.annotated_entry_count`;
`n = input_count` is the number of distinct resolved input accessions, all of which
must belong to that universe. Term membership `K` and input hits `k` come from the
same snapshot. Unannotated/out-of-proteome inputs fail explicitly rather than
being silently removed. Every background term is tested, including zero-hit terms.
The full entry count remains in `universe.proteome_entry_count` for auditing.
`unique_gene_id_count` counts available NCBI GeneID cross-references, not all genes;
`entries_without_gene_id` makes incomplete cross-reference coverage visible.

The shared protein resolver uses exact accession, then exact primary gene symbol,
then exact synonym, with the existing reviewed-entry preference within a symbol
priority. Multiple equally ranked matches remain ambiguous. Inputs are resolved
within the explicit numeric taxon. Repeated identifiers and aliases count once,
while every original input and canonical mapping is retained in the cached result.
Taxon-specific resolver and background cache keys isolate organisms and exclude
old human-only snapshots.

Errors are JSON objects under `error`. Missing context returns `organism_required`;
unresolved/ambiguous identifiers return `resolution_failed` with an `identifiers`
array containing per-identifier reasons. Upstream resolution failures use HTTP 502;
identifier/context failures use HTTP 400. No partial analysis is produced.
Inputs without usable GO BP membership return `outside_annotation_universe`.

See [the background audit](docs/enrichment-background-audit.md) for the original
147,520-entry population, its composition, and the rationale for this universe.

Both input and background term membership come from the **same complete
UniProt `go_p` snapshot**. This reuses the existing UniProt client and avoids
mixing input QuickGO annotations with a different background release or making
one QuickGO request per proteome entry. Existing QuickGO lookup behavior is
unchanged. The analysis uses UniProt's positive BP associations with all evidence,
without adding ancestor propagation; it excludes the BP root `GO:0008150`.
It therefore analyzes the associations supplied by UniProt, and is not an
ancestor-expanded GO analysis. Unannotated inputs are rejected before calculation.

A one-sided hypergeometric survival probability, `P(X >= observed_hits)`, is
computed with `statrs`, equivalent to the enrichment tail of Fisher's exact test.
Benjamini–Hochberg adjustment covers **every non-root BP term observed in the
background**, including terms with zero input hits (p = 1). Membership is a set,
so duplicate annotation records cannot inflate counts.

The first request downloads the paginated background (500 records per page,
at most 1,000 pages), which can take several minutes. Subsequent requests reuse
the complete snapshot. Page totals, unique accession counts and UniProt releases
are checked; HTTP errors, malformed data, truncation or release changes fail with
HTTP 502 and never cache a partial background. No automatic retries are added.
Background loading shares the existing concurrency semaphore. Simultaneous cold
requests can duplicate the load, so warm the cache with one request first.

Snapshots and results use a small `analysis_cache` table in the existing SQLite
file and the existing `CACHE_TTL_DAYS`. They survive restart, expire from write
time, and are not served after expiration. Missing/expired result IDs return
HTTP 404; malformed requests return 400. TTL zero disables retrieval, including
immediate result retrieval. Expired rows have no automatic cleanup, consistent
with the existing cache. No Reactome, KEGG enrichment, GSEA, hierarchy reduction,
network analysis, group comparison or visualization is added by this stage.

References: [UniProt API queries](https://www.uniprot.org/help/api_queries),
[UniProt return fields](https://www.uniprot.org/help/return_fields),
[UniProt proteomes](https://www.uniprot.org/help/proteome), and
[statrs Hypergeometric](https://docs.rs/statrs/0.18.0/statrs/distribution/struct.Hypergeometric.html).
