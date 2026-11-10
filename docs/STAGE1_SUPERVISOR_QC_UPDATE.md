# Stage 1: response to supervisor's coverage-QC feedback

Status: Exploratory technical update, not biological validation.

## Dataset

- ENCODE: ENCSR477GZK; HCT116; GRCh38.
- Region: chr10:20,000,000–21,000,000.
- Resolution: 500 bp; measurement window: 32 kb.
- Scan step: 500 bp; 1,937 windows.
- Contact matrix: ENCFF573OPJ, MAPQ-thresholded.

## Coverage-screening findings

| Region | Original calls | Screened calls | Classification |
|---|---:|---:|---|
| R1 | 23 | 23 | no_coverage_warning |
| R2 | 2 | 2 | no_coverage_warning |
| R3 | 8 | 0 | severe_low_support |
| R4 | 11 | 0 | severe_low_support |
| R5 | 49 | 0 | consecutive_zero_support |

Coverage rules are provisional and were applied consistently
to the complete scan. They are not biologically calibrated.

## Umap-masked, coverage-screened full-scan results

| Umap condition | R1 calls | R2 calls | Other calls | Total |
|---|---:|---:|---:|---:|
| unmasked | 23 | 2 | 0 | 25 |
| k50_50 | 4 | 0 | 0 | 4 |
| k50_80 | 0 | 0 | 0 | 0 |
| k50_95 | 0 | 0 | 7 | 7 |
| k100_50 | 23 | 2 | 0 | 25 |
| k100_80 | 16 | 0 | 0 | 16 |
| k100_95 | 4 | 0 | 0 | 4 |

Counts above include new scan windows within existing
candidate footprints, not only originally called windows.

## Read-length investigation

- Eight upstream FASTQs traced through two BAM files.
- 200 raw reads sampled per FASTQ (1,600 total).
- Observed raw-read median: 113 bp.
- Observed range: 12–241 bp.
- Raw reads have variable lengths.
- Effective mapped read-length distribution remains unknown.
- k50 and k100 remain sensitivity assumptions.

## Interpretation

1. R5 exhibits severe spatial contact-support gaps and is
   withheld by the provisional coverage screen.
2. R3 and R4 also have low-support warnings.
3. R1 and R2 survive the coverage screen but exhibit
   sensitivity to sequence-mappability assumptions.
4. No retained region is considered biologically validated.
5. A zero-contact run does not establish a genomic deletion.

## Remaining work

- Calibrate coverage thresholds beyond the initial simulations.
- Examine effective aligned read lengths where feasible.
- Test an independently justified mappability policy.
- Evaluate independent biological replication.
- Review the candidate classification with Dr. Shamim.
