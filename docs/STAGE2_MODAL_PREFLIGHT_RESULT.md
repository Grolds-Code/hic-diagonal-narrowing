# Stage 2 Modal pooled-data preflight: verified result

**Status:** PASS. Computational/source-provenance verification only;
biological replication remains untested.

- **Run date:** 2026-11-16 (Africa/Nairobi)
- **Modal run:** [ap-ZBeO1IWI6zXvsXZEXhFyhc](https://modal.com/apps/groldotieno97/main/ap-ZBeO1IWI6zXvsXZEXhFyhc)
- **Repository branch:** `research/stage2-replicate-audit`, commit `8261fcf`
- **Script:** `modal_hct116_stage2_replicate_preflight.py`
- **Local JSON:** `research_artifacts/stage2/hct116_stage2_modal_replicate_preflight_v1.json` within the user's **separate local worktree** `~/projects/hic-stage2-modal`

## Results copied from the completed Modal log

| Check | Verified output |
|---|---|
| Execution | Modal app completed successfully; `Stage 2 Modal preflight: PASS` |
| Contact source SHA256 | `fb01bc848b4fecb342a450b8d87165393a89de019fc27c82e72d8519c07f5de5` |
| Stage 1 scan rows independently reproduced by code | 1,937 |
| Provisional retained call windows | 25 |
| R1 footprint | chr10:20,258,500–20,336,000 |
| R1 retained overlapping scan windows | 23; W80 7,000 bp |
| R2 footprint | chr10:20,427,000–20,459,500 |
| R2 retained overlapping scan windows | 2; W80 7,000 bp |
| ENCODE processed products | 12 |
| ENCODE contact-matrix candidates | 4 |
| Biological replicate analysis | NOT TESTED |
| Large sequencing/contact files downloaded | NONE |

The code verifies the existing pooled HCT116 extracted NPZ,
reproduces the archived W80 and contact totals across the full scan, and
rechecks the **provisional** coverage-screen calls. The Modal function
reads the prior Volume without modifying it.

The CSV/JSON output was saved to the user's local detached worktree,
not automatically committed or uploaded to GitHub. The results here
are transcribed from the user's successful Modal console output.
This document does not claim to have independently inspected the
full JSON or any separately processed replicate contact matrices.

## Interpretation

1. Stage 1 provenance and code/data reproducibility: **verified**.
2. Coverage screen: **reproduced**, but thresholds remain provisional.
3. R1/R2 candidate narrowing: **provisional**, with mappability sensitivity.
4. Biological replicate concordance: **not measured**.
5. A single pooled contact extraction cannot stand in for two independent
   biological replicates. Overlapping scan windows are not independent tests.

## Next meaningful Stage 2 step

We must determine whether ENCODE biological-replicate BAM products
`ENCFF624NZV` and `ENCFF583BZD` contain a recoverable Hi-C
contact-pair representation with sufficient provenance to compute
**separate, comparable** R1/R2 contact profiles. Do not treat general
BAM alignment records or template mate positions as validated Hi-C
ligation contacts without verifying the ENCODE pipeline's pairing
semantics, quality filters, duplicate handling and product lineage.

Limit subsequent exploratory probes to metadata / bounded data reads
until suitable extraction methods and compute costs are established.
Only then proceed to separate-replicate W80, coverage, mappability and
pre-specified control comparisons. Preserve all existing local
Stage 2 matrices, plots and manifests.
