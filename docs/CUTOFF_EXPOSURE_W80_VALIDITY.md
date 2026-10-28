# Statistical validity of cutoff-specific W80 confidence intervals

## 1. Population estimand

For contact separations d = 1, ..., D, define:

- K_d > 0: fixed eligible genomic-pair count.
- mu_d >= 0: population mean contact intensity per eligible pair.
- S = sum(mu_d) > 0.

The population cumulative intensity profile is

    F(d) = sum(mu_j, j <= d) / S.

The population W80 is

    w = min{d : F(d) >= 0.80}.

This is a contact-distance statistic, not a direct measurement
of three-dimensional chromatin compaction.

## 2. Sampling assumptions

Assume independent band totals:

    Y_d ~ Poisson(K_d * mu_d).

Let N = sum(Y_d).

Conditional on N = n > 0, the band counts are multinomial,
with probabilities proportional to K_d * mu_d.

Consequently, the observed cumulative count C_d satisfies

    C_d | N=n ~ Binomial(n, Q_d),

where

    Q_d = sum(K_j * mu_j, j <= d)
          / sum(K_j * mu_j, all j).

No equal-exposure assumption is necessary.

## 3. Cutoff-specific exposure bounds

For every cutoff d < D, partition the distances into:

- Prefix: separations <= d.
- Tail: separations > d.

Let the minimum and maximum exposures within each
partition be K_prefix_min, K_prefix_max,
K_tail_min, and K_tail_max.

Set t = 0.80.

If F(d) >= t, then

    Q_d >= p_low(d),

where

    p_low(d) =
        t * K_prefix_min /
        (t * K_prefix_min + (1-t) * K_tail_max).

If F(d) <= t, then

    Q_d <= p_high(d),

where

    p_high(d) =
        t * K_prefix_max /
        (t * K_prefix_max + (1-t) * K_tail_min).

These inequalities follow from bounding the
exposure-weighted prefix and tail intensities.

Neither probability must straddle 0.80.

## 4. Confidence interval construction

For confidence level 1-alpha and observed N=n>0:

    a = alpha/2

At each cutoff d < D, calculate

    r_low(d) = BinomialPPF(a; n, p_low(d))
    r_high(d) = BinomialISF(a; n, p_high(d)).

Define L as the earliest cutoff satisfying

    C_d >= r_low(d).

If none exists, set L = D.

Define U as one plus the latest cutoff satisfying

    C_d <= r_high(d).

If none exists, set U = 1.

The resulting interval is [L, U].

The construction does not assume that the
cutoff-specific thresholds are monotonic.

## 5. Conditional coverage argument

Let w be the true population W80.

When w < D, F(w) >= t implies

    Q_w >= p_low(w).

Binomial stochastic ordering and the lower-tail
quantile definition imply

    P(C_w < r_low(w) | N=n) <= alpha/2.

If C_w >= r_low(w), the cutoff w is an admissible
lower-bound candidate. Therefore L <= w.

Thus

    P(L > w | N=n) <= alpha/2.

When w > 1, F(w-1) < t implies

    Q_(w-1) <= p_high(w-1).

The corresponding upper-tail bound implies

    P(C_(w-1) > r_high(w-1) | N=n) <= alpha/2.

If C_(w-1) <= r_high(w-1), the cutoff w-1 is
an admissible upper-bound candidate. Therefore U >= w.

Thus

    P(U < w | N=n) <= alpha/2.

At the boundary values w=1 and w=D, the
respective one-sided error is impossible.

Applying the union bound yields

    P(L <= w <= U | N=n) >= 1-alpha

for every n>0 under the specified sampling model.

Only the two cutoffs adjoining the true W80 are
needed for this argument, not simultaneous
testing across every separation.

## 6. Limitations

This guarantee depends on:

- Independent Poisson contact-band totals.
- Fixed, positive eligible-pair exposures.
- Raw nonnegative integer contact counts.
- A positive observed total contact count.
- A population W80 defined by mean intensity.

The implementation abstains when no contacts are observed.

The guarantee does not automatically transfer to:

- Biological Hi-C with contact dependence or overdispersion.
- Data-dependent missing-bin masks.
- Normalized or transformed contact matrices.
- Selected genomic candidates or multiple-window inference.
- Comparisons of W80 between two genomic windows.

Independent biological validation remains necessary.

## 7. Computational evidence

The reusable implementation matched the independent
experimental implementation on 7,200 synthetic trials,
with zero mismatches.

A separate exact-enumeration audit evaluated 144
small-model configurations. All achieved at least
95% conditional coverage.

These tests support correctness but do not replace
the statistical argument or validate its assumptions
for biological Hi-C.

## 8. Implementation and verification

Implementation:
    src/hic_diagonal_narrowing/poisson_uncertainty.py

Function:
    cutoff_exposure_w80_interval_from_bands

Tests:
    tests/test_cutoff_exposure.py

Experiments:
    experiments/compare_cutoff_exposure_rank.py
    experiments/verify_cutoff_exposure_rank_exact.py
    experiments/verify_cutoff_exposure_implementation.py
