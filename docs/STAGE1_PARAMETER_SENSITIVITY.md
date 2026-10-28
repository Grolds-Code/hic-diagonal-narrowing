# Stage 1: HCT116 parameter-sensitivity analysis

Status: Exploratory technical validation.

## 1. Objective

Assess how the original HCT116 narrowing calls depend on
the provisional width-estimation and scanning parameters.

This is a sensitivity analysis, not a biological hypothesis
test or a post-selection significance calculation.

## 2. Experimental design

- Dataset: HCT116 chr10:20,000,000–21,000,000.
- Resolution: 500 bp.
- Window sizes: 16, 32, and 48 kb.
- Maximum contact separations: 16, 24, and 32 kb.
- Cumulative-intensity fractions: $q \in \{0.70, 0.80, 0.90\}$.
- Genomic scan step: 8 kb.
- Original narrowing-score threshold: $S(W) > 0.2$.
- Maximum separation must not exceed window length.
- Valid parameter configurations: 21.
- Original candidate centers were held fixed.

All profiles use midpoint-assigned contacts and the
original detector's non-overlapping reference rule.

Changing the cumulative-intensity fraction changes the
width estimand. Results across $W_{70}$, $W_{80}$, and $W_{90}$ therefore
represent sensitivity to the estimand definition.

## 3. Baseline verification

The original 32-kb window, 32-kb maximum separation,
and 80% cumulative fraction reproduced:

- 122 genomic windows.
- Five called windows.
- Two merged candidate regions.
- Region A: $W_{80} = 14$ bins.
- Region B: $W_{80} = 13$ bins.

Baseline reconstruction: PASS.

## 4. Candidate-center sensitivity

| Window size | Configurations | Region A calls | Region B calls |
|---|---:|---:|---:|
| 16 kb | 3 | 1 | 2 |
| 32 kb | 9 | 3 | 8 |
| 48 kb | 9 | 0 | 1 |
| **Total** | **21** | **4** | **11** |

Region B was detected in eight of nine configurations
using the original 32-kb window length.

At a 48-kb window length, Region A was not detected
at its original center in any tested configuration,
and Region B was detected in one of nine.

These are descriptive configuration counts, not
independent-trial probabilities.

## 5. Interpretation

Region B exhibits greater parameter stability than
Region A within this specific sensitivity grid.

Neither region is detected consistently across all
tested parameter combinations.

The results show that window size, maximum contact
separation, and cumulative-intensity definition can
materially affect candidate detection.

A candidate not called at its original center under
a different window size may shift to an adjacent window.
This analysis does not establish the absence of nearby
narrowing signals.

The selected parameter grid is exploratory and is not
a biologically calibrated uncertainty distribution.

## 6. Limitations

- Candidate centers were selected from the baseline data.
- Reference medians were recomputed for each configuration.
- Different width fractions and separation cutoffs define
  different estimands.
- Configurations are statistically dependent.
- Original genomic QC remains unassessed.
- Independent missing-bin validity is not established.
- Biological replicate confirmation is unavailable.
- No claim of biological significance is supported.

## 7. Stage 1 decision

Report both regions as exploratory narrowing candidates.

Region B has stronger technical and parameter stability,
but neither region can be described as universally robust
or biologically validated.

The original 32-kb / 32-kb / 80% configuration remains
the documented baseline, not an empirically established
optimal configuration.

## 8. Reproducibility

- Experiment:
  `modal_hct116_parameter_sensitivity.py`
- Results:
  `research_artifacts/stage1/hct116_parameter_sensitivity_v1.json`
- Source SHA-256:
  `fb01bc848b4fecb342a450b8d87165393a89de019fc27c82e72d8519c07f5de5`
- Baseline profile SHA-256:
  `b9205d00020568e038be3bbb42dfde3644d00cb86dcff127c1ea7f0abedf662f`
