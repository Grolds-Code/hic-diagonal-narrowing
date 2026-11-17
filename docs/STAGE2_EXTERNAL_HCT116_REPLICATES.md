# Stage 2: alternate independent HCT116 validation dataset

**Discovery date:** 2026-11-17.  
**Status:** Metadata/source discovery only; no downloaded matrices, no W80 concordance.

## Why a new dataset

The original ENCODE experiment `ENCSR477GZK` has two separately labeled
biological-replicate BAMs, but the two `.hic` products with explicit
`derived_from` both use those two BAMs. Two older submitted `.hic`
products were initially unlabeled; their submitted paths include
`untreated/ENCODE_18cell/combined/total/inter.hic` and
`untreated/ENCODE_18cell/combined/total/inter_30.hic`, strongly indicating
combined contact maps. The file names correspond to different MAPQ
thresholds, **not independent replicates**.

ENCODE explicitly describes a `mapping quality thresholded contact
matrix` as a contact matrix filtered to MAPQ >= 30:
https://www.encodeproject.org/hic/

## New candidate: GEO GSE310051

The series **GSE310051**, *Hyperosmotic stress induces a large-scale
rewiring of 3D chromatin interactions [Hi-C]*, has four individual
**HCT116 wild-type, untreated control** samples at the 0h timepoint.
Separate replicate `inter_30.hic` files are listed by NCBI GEO,
assembled to **hg38/GRCh38**, and no sorbitol treatment was applied
to these four controls.

| Biological sample | GEO sample | Individual matrix | GEO-listed approximate size |
|---|---|---|---|
| Control 1 | `GSM9287304` | `GSM9287304_STRS_HCT116_WT_control_0h_1_1_inter_30.hic` | 2.1 GB |
| Control 2 | `GSM9287305` | `GSM9287305_STRS_HCT116_WT_control_0h_2_1_inter_30.hic` | 2.1 GB |
| Control 3 | `GSM9287306` | `GSM9287306_STRS_HCT116_WT_control_0h_3_1_inter_30.hic` | 2.1 GB |
| Control 4 | `GSM9287307` | `GSM9287307_STRS_HCT116_WT_control_0h_4_1_inter_30.hic` | 2.3 GB |

Primary sample sources:
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE310051
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9287304
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9287305
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9287306
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9287307

**Experimental caveat:** These controls used **in situ MboI
restriction enzyme Hi-C**, while ENCODE `ENCSR477GZK` specifies
**MNase-based intact Hi-C**. This makes them *independent biological
experimental data*, but the datasets cannot be treated as identical
technical measurements. Differences in restriction/proximity-ligation
protocol, depth, mappability, normalization and read filtering could
change W80 at 500 bp even for unchanged biology.

### Mandatory feasibility checks before analysis

1. Check URLs, range-request capability, `.hic` format version and
   accessible resolution metadata **without downloading whole files**.
2. Verify that 500 bp really exists in *each* map and evaluate
   chr10:20–21 Mb local contact occupancy and marginal read support.
   Presence of a 500-bp zoom level alone does not imply adequate power.
3. Validate the extracted per-distance contact representation against
   the original Stage 1 *midpoint-assigned*, 32-kb-window, 32-kb maximum
   separation estimator. Use raw `observed` contact counts with no
   arbitrary balanced-vs-unbalanced substitution.
4. Pre-register the R1/R2 footprints, all negative-control windows,
   mappability and missing-bin QC, and depth-matching sensitivity
   **before viewing candidate results in the validation dataset**.
5. Treat four control samples as four distinct biological libraries
   from the *new* study, not as replicate-specific reconstructions
   from the original ENCODE experiment.
6. Report any cross-protocol failure to reproduce. Do not retrofit
   thresholds or move candidate regions to fit external datasets.

## Decision

**Prefer this independent, lower-storage validation route over
immediate reprocessing of ~942 GB of ENCODE replicate BAMs.** The
decision is conditional on 500-bp support and adequate local depth.
No biological replication claim can be made yet.
