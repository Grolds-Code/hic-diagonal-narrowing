# Stage 2 processed Hi-C lineage audit

**Date:** 2026-11-17 (Africa/Nairobi)  
**Status:** Metadata-only audit PASS. No biological replicate W80 test performed.

**Source:** User's Step 401 local terminal output from the script
`experiments/stage2_processed_lineage.py` at commit `285166c`.
The detailed local JSON is at
`~/projects/hic-stage2-modal/research_artifacts/stage2/hct116_stage2_processed_lineage_v1.json`.
The local JSON was not committed or independently reviewed here.

## Verified file-level metadata

| File | Format and output type | Size (bytes) | Biological replicate fields | Direct ENCODE `derived_from` |
|---|---|---:|---|---|
| `ENCFF573OPJ` | `hic`, MAPQ-thresholded contact matrix | 43,687,472,606 | 1, 2 | `ENCFF583BZD`, `ENCFF624NZV` |
| `ENCFF726YRX` | `hic`, contact matrix | 46,197,948,255 | 1, 2 | `ENCFF583BZD`, `ENCFF624NZV` |
| `ENCFF743HCX` | `hic`, MAPQ-thresholded contact matrix | 107,243,001,887 | none provided | none provided |
| `ENCFF876NPA` | `hic`, contact matrix | 110,318,183,165 | none provided | none provided |
| `ENCFF498SJT` | `pairs`, pairs | 13,736,762,544 | 1, 2 | `ENCFF583BZD`, `ENCFF624NZV` |

`ENCFF573OPJ` and `ENCFF726YRX` share the same step-run reference
`/analysis-step-runs/ce24e902-47e8-478e-ae18-453925eb7477/`.
The `.pairs` file uses a different processing step-run:
`/analysis-step-runs/ed197695-0ed6-4e98-8d08-f44baf451858/`.
No `analysis`, `pipeline` or `step_run` reference was supplied
for the two lineage-unresolved matrices in the inspected fields.

**Derived conclusion:** The three files with explicit ancestry
(two `.hic` files and one `.pairs` file) are pooled-source products,
not evidence of two independently processed biological-replicate maps.

**Non-conclusion:** Missing replicate labels and missing direct parents
in `ENCFF743HCX` and `ENCFF876NPA` do **not** establish whether
those matrices are pooled or replicate-specific. Their larger sizes
are not evidence of independence. Neither should be downloaded solely
to decide provenance.

## Scientific status

- The pooled Modal 500-bp W80 preflight passed all 1,937 archived scan
  measurements, but that is a **code/data regression**, not biological replication.
- R1 has 23, R2 has 2 provisionally retained overlapping calls;
  both remain mappability-sensitive.
- No separate biological-replicate contact matrix has yet been
  identified with sufficient lineage.
- Stage 2 is **not complete**.

## Recommended immediate step

Inspect the **already downloaded small JSON metadata** for the two
lineage-unresolved matrices, especially `assembly`, `dataset_reference`,
`submitted_file_name`, `accession_labels`, `status`, and any
available processing identifiers. If these are inconclusive, inspect the
ENCODE experiment/analysis-step or associated GEO provenance, and consider
contact-pair outputs from a separately documented experiment.

A full reprocessing of replicate BAMs is **not** authorized merely by
this metadata audit. Earlier verified BAM file sizes were 159,057,417,163
and 782,504,438,499 bytes (approximately 942 GB total). Before any
large Modal transfer or computation, require a documented pair-reconstruction
protocol, quality criteria, resource budget and validation plan.

## Step 402: submitted-file paths clarify previously unresolved matrices

The user inspected the locally saved metadata JSON. These **released,
GRCh38** files in `ENCSR477GZK` have the following submitted names:

- `ENCFF743HCX`: `.../HCT-116/untreated/ENCODE_18cell/combined/total/inter_30.hic`,
  alias `erez-aiden:intact-HCT116-hic30map`
- `ENCFF876NPA`: `.../HCT-116/untreated/ENCODE_18cell/combined/total/inter.hic`,
  alias `erez-aiden:intact-HCT116-hic1map`

`combined/total` is strong **file-path evidence** for pooled rather
than replicate-separated maps, but no corresponding direct
`derived_from` ancestry was supplied for these two files. The Juicer
filenames `inter.hic` versus `inter_30.hic` conventionally refer to
MAPQ >0 and MAPQ >=30 filtering; **these are not two biological
replicates**.

Primary reference:
https://www.encodeproject.org/hic/
Juicer documentation:
https://github.com/aidenlab/JuicerTools/blob/master/Juicer.README

Therefore the four analyzed `.hic` products do not supply a
verified **pair of independent replicate-specific contact maps** for
the original ENCODE sample. Before reprocessing nearly a terabyte of
replicate alignments, investigate smaller externally published
HCT116 untreated control replicates, as tracked in
[STAGE2_EXTERNAL_HCT116_REPLICATES.md](STAGE2_EXTERNAL_HCT116_REPLICATES.md).
