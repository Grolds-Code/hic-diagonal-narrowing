# Stage 1 QC decision record

Status: Exploratory technical analysis; formal genomic QC not passed or failed.

## 1. Dataset provenance

- Source: HCT116 chr10:20,000,000–21,000,000, 500-bp resolution.
- Source SHA-256: `fb01bc848b4fecb342a450b8d87165393a89de019fc27c82e72d8519c07f5de5`.
- Profile SHA-256: `b9205d00020568e038be3bbb42dfde3644d00cb86dcff127c1ea7f0abedf662f`.
- Contact normalization: `NONE`.
- Contact-pair ownership: midpoint-based.
- Number of evaluated windows: 122.
- Recorded QC status: `unassessed` for every window.
- Independent missing-bin status: `not_assessed`.

The input checksums establish data consistency, not biological
quality or validity of a missing-bin mask.

## 2. Descriptive measurement support

Across the evaluated windows:

- Contact totals: 1737 to 8102.
- Pair occupancy: 0.1765 to 0.5300.
- W80: 13 to 20 bins.
- All W80 values are finite and positive.
- All windows contain positive recorded contact signal.

These observations demonstrate computable descriptive profiles.
They do not constitute a formal QC pass.

## 3. Exploratory candidate measurements

| Quantity | Region A | Region B |
|---|---:|---:|
| Peak start (bp) | 20296000 | 20736000 |
| W80 (bins) | 14 | 13 |
| W80 (bp) | 7000 | 6500 |
| Total contacts | 2478 | 4695 |
| Pair occupancy | 0.235840 | 0.303711 |
| Supported distance bands | 64 | 64 |
| QC status | unassessed | unassessed |

A supported distance band means that signal was observed for
that band, not that all genomic bins are independently verified.

## 4. What has been validated

- Original W80 profiles independently reproduced from source contacts.
- Midpoint-based contact geometry audited.
- Synthetic coverage-loss and depth-stress experiments performed.
- Real HCT116 technical thinning and joint candidate-score sensitivity performed.
- Model-conditional uncertainty studied analytically and computationally.
- Spatial and smooth exploratory null simulations performed.
- Genomic mappability sensitivity investigated.

These validations address specified technical failure modes.
They do not establish biological significance.

## 5. What is not established

1. An independently verified missing-bin or valid-bin mask for this Hi-C dataset.
2. A biologically calibrated minimum contact-count or occupancy threshold.
3. Formal QC acceptance criteria with demonstrated error-control properties.
4. Independent biological replicate stability of these specific candidate scores.
5. Biological significance after selection across the genomic scan.
6. Precise structural boundaries or direct chromatin-compaction measurements.

In particular, a narrow W80 alone cannot be interpreted as
physical chromatin compaction.

## 6. Stage 1 reporting policy

- Retain the recorded `unassessed` status for all 122 windows.
- Report measured support separately from model-based uncertainty.
- Label Regions A and B as exploratory candidate narrowing signals.
- Report Region B's stronger technical persistence without claiming
  independently validated biological superiority.
- Identify the original detector's QC acceptance policy as provisional.
- Do not introduce retrospective acceptance thresholds.
- Do not report exploratory Monte Carlo null fractions as
  biologically calibrated p-values.
- Treat missing-bin validity and biological replication as
  explicit limitations and future validation requirements.

## 7. Decision

The current evidence is sufficient to report a reproducible
methodological case study with exploratory candidate detection.

It is not sufficient to promote individual windows to
biologically validated narrowing calls.

Formal QC implementation and acceptance calibration remain
outstanding for confirmatory or genome-wide deployment.
