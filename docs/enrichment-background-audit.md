# Stage 1 resolver and GO BP universe audit

Audited 2026-09-28 (America/Chicago), UniProt release **2026_03**, released
2026-09-02. Counts below come from a completed live download, not an estimate.

## Root causes and changes

The shared resolver searched `gene_exact` and treated a primary gene name and a
synonym as equal candidates. In taxon 9606, `ATR` matched three reviewed entries:
Q13535 (primary gene ATR), P20848 (SERPINA2, synonym ATR), and Q9H6X2 (ANTXR1,
synonym ATR). Its ambiguity error was correct under that old ranking policy.
The tool failed the batch; the later manual translation and omission happened
outside the service. Enrichment also hardcoded the name Homo sapiens and exposed
no organism context.

The shared resolver now accepts numeric taxon context, validates returned taxon
IDs, and ranks exact primary symbols ahead of exact synonyms. The existing
reviewed-entry preference applies within that priority. Equal-priority ambiguity
still fails. Enrichment aggregates errors before running any analysis and retains
every original input-to-canonical mapping. The regression uses exactly TP53, ATM,
ATR, BRCA1, BRCA2, CHEK1, CHEK2, RAD51, MDM2, CDKN1A.

The old background loader requested `accession,go_p` from
`proteome:UP000005640 AND organism_id:9606`. It inserted each accession into a set
**before** checking annotations, verified its size against `x-total-results`, and
passed that entire set's length to `Hypergeometric::new`. Therefore the old
population parameter was **N = 147,520**, including unannotated records. It did
not count GO rows, isoform FASTA sequences, input proteins, or unique genes.

The smallest statistical change retains the accession unit and conditions the
universe on having usable BP membership. No gene collapsing, reviewed-only
filter, or fixed approximately 20,000-gene population was introduced. The
reference proteome is now discovered from the requested exact taxon. Missing,
multiple or incompletely enumerated reference proteomes fail explicitly.

## Composition of the original population

| Measurement | Count |
|---|---:|
| Rows / unique primary UniProtKB accessions | 147,520 |
| Reviewed entries | 20,416 |
| Unreviewed entries | 127,104 |
| Entries with non-root GO BP associations | 36,939 |
| Annotated reviewed entries | 17,164 |
| Annotated unreviewed entries | 19,775 |
| Entries without usable GO BP associations | 110,581 |
| Distinct NCBI GeneID cross-references | 19,416 |
| Entries with at least one GeneID | 31,870 |
| Entries without a GeneID | 115,650 |
| GeneIDs associated with multiple entries | 6,161 |
| Distinct primary gene-name tokens | 20,699 |
| Entries without primary gene names | 952 |
| Expanded `accession-N` isoform rows | 0 |

All 147,520 accession strings are unique. Each represents a canonical **entry**;
this does not mean one canonical protein per gene. Multiple entries represent
products of the same gene: this download contains 18 entries named TP53, 8 named
ATR, and 12 named PTEN, with one reviewed entry in each group. The gene-ID
coverage is incomplete, so 19,416 is the number of represented, mapped GeneIDs,
not an exact census of all genes. Similarly, 20,699 distinct gene-name tokens
must not be represented as an authoritative unique-gene count.

The TSV query does not expand curated isoforms such as `P04637-2` into extra rows.
Unreviewed entries can nevertheless represent predicted alternative products of
already represented genes. UniProt describes these distinctions in its
[proteome documentation](https://www.uniprot.org/help/proteome) and its primary
[guide to the human proteome](https://pmc.ncbi.nlm.nih.gov/articles/PMC4761109/).
Those sources explain the data model; the counts above are from the live 2026_03
download, not the historical paper. This audit does not claim that repeated gene
labels imply identical sequences, nor does it attempt sequence-level deduplication.

## Statistical population after the fix

For `organism_taxon=9606`, discovery returns UP000005640. For any supported taxon:

1. Discover the unique reference proteome for the exact taxonomy ID.
2. Download its entries under `proteome:<id> AND organism_id:<taxon>`.
3. Build positive GO BP membership from the same complete `go_p` snapshot,
   excluding GO:0008150 and performing no additional ancestor propagation.
4. Define eligible accessions as the union of those term-membership sets.

For this release, **N = 36,939**. Each distinct resolved input accession must be
in that universe. An unresolved, ambiguous, wrong-taxon, out-of-proteome, or
unannotated input causes an explicit error; none is silently removed. Successful
analysis of the ten-symbol regression has **n = 10**. For each term, K is its
eligible background membership and k is its input membership. The calculation
uses P(X >= k), followed by BH across all non-root background terms, including
zero-hit terms.

This is annotation-conditioned **protein-entry** enrichment. It asks about term
representation among GO BP-annotated canonical entries in the reference proteome.
It is not gene-level enrichment and does not claim entries from the same gene
are independent genes. This distinction remains relevant for interpreting lists
entered as gene symbols. A gene-level analysis would require a separate,
explicit, complete gene mapping and gene-level membership policy; changing N to
20,000 while leaving accession-level K and k would be inconsistent.

The cached full result records the taxon, reference proteome ID, background query,
UniProt release, annotation/selection policies, analysis unit, full and annotated
entry counts, review counts, GeneID coverage, exact `eligible_accessions`, original
inputs, and canonical mappings. The snapshot cache is versioned and taxon-scoped;
old human-only snapshots cannot supply the new universe.

## Reproducing and inspecting the audit

The completed audit files in this workspace session are:

- `/tmp/protein-background-audit.tsv`: full original population, with review,
  primary gene name, GeneID and GO BP fields.
- `/tmp/protein-background-audit.headers`: release headers.
- `/tmp/protein-background-audit-counts.json`: measured counts.
- `/tmp/protein-background-eligible-accessions.txt`: every one of the 36,939
  entries with at least one non-root GO BP association, one accession per line.

The TSV SHA-256 is
`19612c6d01e73165393fa2a94550503c9e578fe34569f5c706ef415fba080ddb`.

```sh
curl -sS --max-time 900 --compressed -G \
  'https://rest.uniprot.org/uniprotkb/stream' \
  --data-urlencode 'query=proteome:UP000005640 AND organism_id:9606' \
  --data-urlencode 'format=tsv' \
  --data-urlencode 'fields=accession,reviewed,gene_primary,xref_geneid,go_p' \
  -D /tmp/protein-background-audit.headers \
  -o /tmp/protein-background-audit.tsv
```

Count distinct `Entry` values, split `GeneID` on semicolons and primary gene names
on whitespace. An entry is eligible if at least one semicolon-separated BP term
has an ID other than GO:0008150. Counts and checksums can change with a later
UniProt release. Production uses paginated search with total/release checks,
rather than the audit's streaming download.

The offline integration regression includes the actual ATR alias collision,
reviewed/unreviewed candidates, all ten inputs, provenance, duplicate input
handling, taxon isolation, aggregated failures and annotation eligibility. The
opt-in live integration runs the actual shared resolver, reference-proteome
discovery, full background download and statistical calculation:

```sh
cargo test
cargo test --test enrichment live_uniprot_exact_ten_protein_set -- --ignored --nocapture
cargo clippy --all-targets --all-features -- -D warnings
```

## Validation results

- Full suite: **38 passed**, zero failures; the external-network test is opt-in.
- Live integration, run separately: **passed**, all ten identifiers included;
  release 2026_03, input_count=10, N=36,939, proteome_entry_count=147,520.
  The complete cold-background run took 495.98 seconds.
- `cargo clippy --all-targets --all-features -- -D warnings`: passed.
- `cargo fmt --check` and `git diff --check`: passed.
