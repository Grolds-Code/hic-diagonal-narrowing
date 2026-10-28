# Stage 1: Hi-C Diagonal Narrowing Detection

Status: Provisional research specification (v0.1)

## 1. Scientific objective

Identify genomic regions where the local Hi-C contact
diagonal is substantially narrower than an appropriate
reference, while distinguishing biological structure
from coverage limitations and technical artifacts.

## 2. Initial mathematical approach

For a genomic window $W$:

1. Select valid cis-contact pairs whose genomic midpoint
   falls inside $W$.
2. Restrict contact separation to a fixed range:
   $0 < d \le D_{\max}$.
3. Calculate mean observed contact intensity at each
   separation $d$, accounting for the number of valid
   bin pairs available at that separation.
4. Normalize these intensities into a local
   contact-distance distribution.
5. Define width $W_{80}$ as the separation containing 80%
   of the cumulative contact-distance signal.

The 80% threshold, window size, and maximum separation
are provisional parameters requiring sensitivity analysis.

## 3. Candidate narrowing score

Initial proposal:

$$
S(W) = \log_2\left(\frac{W_{80}^{\mathrm{ref}}(W)}{W_{80}(W)}\right)
$$

Positive scores indicate narrower local contact profiles.

The reference may be estimated from neighboring or matched
genomic windows. Reference selection requires validation.

## 4. Essential quality controls

- Minimum informative contact count.
- Missing-bin and invalid-bin masks.
- Distance-dependent contact occupancy.
- Window-edge and chromosome-boundary effects.
- Sensitivity to sequencing depth and downsampling.
- Comparison of raw and appropriately normalized signals.
- Avoid interpreting coverage loss as biological narrowing.

## 5. Initial validation requirements

Before applying the detector genome-wide:

- Recover known narrowing in synthetic contact matrices.
- Reject uniform intensity reductions as narrowing.
- Identify unreliable low-coverage regions.
- Test robustness to missing contacts.
- Examine parameter sensitivity.
- Compare against matched negative controls.

## 6. Stage 1 outputs

- Local diagonal-width profiles.
- Candidate genomic intervals.
- Narrowing scores.
- Coverage and quality-control metrics.
- Reproducible figures.
- Documented Python analysis pipeline.

## 7. Scientific limitations

$W_{80}$ is an operational width statistic, not a direct
measurement of physical chromatin compaction.

Biological interpretation requires additional validation.
No biological significance is assumed from detection alone.
