# Stage 2: ENCODE replicate-data feasibility findings

**Audit date:** 2026-11-16 (UTC)  
**Experiment:** ENCODE `ENCSR477GZK` (HCT116, GRCh38)  
**Method:** Metadata-only automated GitHub Actions audit; no contact or sequencing
data downloaded  
**Evidence:** [GitHub Actions run 43851527908](https://github.com/Grolds-Code/hic-diagonal-narrowing/actions/runs/43851527908),
[JSON artifact](https://github.com/Grolds-Code/hic-diagonal-narrowing/actions/runs/43851527908/artifacts/14218350617)

## Inventory

The ENCODE search returned **28 files**, including **four `.hic` files**,
**two BAM alignment files**, **one `.pairs` file**, **eight FASTQs**,
five `.bedpe` files, three BED files and five bigWigs. The audit
fetched detailed metadata for 12 processed products.

| ENCODE accession | Type / output | Replicate labels | Direct lineage | Interpretation |
|---|---|---|---|---|
| `ENCFF583BZD` | BAM alignments | Biological 2; technical 2_1, 2_2 | Four FASTQs plus index | Available biological-replicate-2 alignments |
| `ENCFF624NZV` | BAM alignments | Biological 1; technical 1_1, 1_2 | Four FASTQs plus index | Available biological-replicate-1 alignments |
| `ENCFF573OPJ` | MAPQ-thresholded `.hic` | Biological 1 **and** 2 | Both BAMs above | Pooled-source contact matrix used in Stage 1 |
| `ENCFF726YRX` | Contact matrix `.hic` | Biological 1 **and** 2 | Both BAMs above | Another pooled-source contact matrix |
| `ENCFF743HCX` | MAPQ-thresholded `.hic` | No replicate label | No direct parents listed | Replicate provenance unresolved; do not assign |
| `ENCFF876NPA` | Contact matrix `.hic` | No replicate label | No direct parents listed | Replicate provenance unresolved; do not assign |
| `ENCFF498SJT` | `.pairs` file | Biological 1 **and** 2 | Both BAMs above | Pooled-source pairs product |

The experiment metadata lists two biological-replicate groups with two
technical replicates each. A biological-replicate-level BAM is not by
itself an independently processed replicate-level Hi-C matrix.

## Research decision

**Replicate-specific alignment data are available, but the inspected
annotated contact matrices do not establish two independently processed
replicate-level matrices.**

Before claiming replicate concordance:

1. Verify the exact BAM alignment schema, contact-pair semantics,
   normalization and chromosome naming used by the ENCODE Hi-C pipeline.
2. Determine whether one can create compatible replicate-specific
   500-bp contact matrices from the two BAMs, ideally preserving the
   Stage 1 midpoint assignment and maximum 32-kb genomic separation.
3. Verify sufficient read depth, bin support and sequence-mappability
   support *per replicate*. The pooled matrix's quality does not transfer
   automatically to either replicate.
4. Compute W80 for R1 and R2 and prospectively specified control windows
   separately for both biological replicate groups. Avoid changing
   acceptance thresholds after looking at the results.
5. Preserve Stage 1's `formal_qc_status=unassessed` until justified
   acceptance criteria are established; report all negative findings.

The two unlabeled `.hic` files should be investigated before ruling out
alternate independent data products; missing lineage in metadata is
not proof that they are pooled or replicate-specific.

## Boundaries

- This is **a feasibility inventory**, not a replication result.
- No contact matrix has been recomputed on this branch.
- The Stage 2 exploratory aggregate and mappability scripts remain
  untracked on the user's WSL machine, outside GitHub.
- No Stage 1 result has been modified.
- Scientific assumptions about deletion, narrowed structure or biological
  significance remain unverified.
