"""Compare global and cutoff-specific exposure-bounded W80 intervals.

Assumptions:
    Independent Poisson contact counts.
    Fixed positive exposure per separation.
    Positive observed total contact count.

No biological Hi-C confidence-interval validation is implied.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binom

D_MAX = 64
DISTANCES = np.arange(1, D_MAX + 1)
N_REPLICATES = 200
ALPHA = 0.05
TARGET = 0.80

mask_rng = np.random.default_rng(20261028)

valid_bins = np.ones(256, dtype=bool)
valid_bins[
    mask_rng.choice(256, size=58, replace=False)
] = False

masked = np.zeros(D_MAX, dtype=np.int64)

for d in DISTANCES:
    rows = np.arange(256 - d)
    midpoint = rows + d / 2

    selected = (
        (midpoint >= 96)
        & (midpoint < 160)
        & valid_bins[rows]
        & valid_bins[rows + d]
    )

    masked[d - 1] = selected.sum()

exposures = {
    "equal": np.full(D_MAX, 40, dtype=np.int64),
    "masked": masked,
    "decreasing": np.rint(
        np.linspace(70, 10, D_MAX)
    ).astype(np.int64),
    "increasing": np.rint(
        np.linspace(10, 70, D_MAX)
    ).astype(np.int64),
}

baseline = 90 * np.exp(-DISTANCES / 18)
narrowing = 90 * np.exp(-DISTANCES / 5)

intensities = {
    "baseline": baseline,
    "narrowing": narrowing,
    "coverage_drop": 0.25 * baseline,
}

truths = {
    "baseline": 28,
    "narrowing": 9,
    "coverage_drop": 28,
}


def global_bounds(counts, exposure):
    n = int(counts.sum())

    if n == 0:
        return None

    kmin = exposure.min()
    kmax = exposure.max()

    p_low = TARGET * kmin / (
        TARGET * kmin + (1 - TARGET) * kmax
    )

    p_high = TARGET * kmax / (
        TARGET * kmax + (1 - TARGET) * kmin
    )

    lower_rank = int(
        binom.ppf(ALPHA / 2, n, p_low)
    )

    upper_rank = int(
        binom.isf(ALPHA / 2, n, p_high)
    ) + 1

    cumulative = np.cumsum(counts)

    lower = (
        1 if lower_rank < 1
        else int(np.searchsorted(cumulative, lower_rank) + 1)
    )

    upper = (
        D_MAX if upper_rank > n
        else int(np.searchsorted(cumulative, upper_rank) + 1)
    )

    return lower, upper


def cutoff_specific_bounds(counts, exposure):
    n = int(counts.sum())

    if n == 0:
        return None

    # For each cutoff d = 1,...,D_MAX-1:
    # prefix = separations <= d
    # tail = separations > d

    prefix_min = np.minimum.accumulate(exposure)[:-1]
    prefix_max = np.maximum.accumulate(exposure)[:-1]

    suffix_min = np.minimum.accumulate(
        exposure[::-1]
    )[::-1][1:]

    suffix_max = np.maximum.accumulate(
        exposure[::-1]
    )[::-1][1:]

    p_low = TARGET * prefix_min / (
        TARGET * prefix_min
        + (1 - TARGET) * suffix_max
    )

    p_high = TARGET * prefix_max / (
        TARGET * prefix_max
        + (1 - TARGET) * suffix_min
    )

    # Exposure-weighted contact probabilities need not
    # straddle TARGET when exposures differ by separation.
    # The required property is ordered, valid probability bounds.
    assert np.all(np.isfinite(p_low))
    assert np.all(np.isfinite(p_high))
    assert np.all(p_low > 0)
    assert np.all(p_high < 1)
    assert np.all(p_low <= p_high)

    lower_thresholds = binom.ppf(
        ALPHA / 2, n, p_low
    )

    upper_thresholds = binom.isf(
        ALPHA / 2, n, p_high
    )

    cumulative = np.cumsum(counts)[:-1]

    # Earliest cutoff compatible with reaching 80%.
    lower_candidates = np.flatnonzero(
        cumulative >= lower_thresholds
    )

    lower = (
        int(lower_candidates[0]) + 1
        if lower_candidates.size
        else D_MAX
    )

    # Latest cutoff where an 80% crossing is not yet
    # established. Using the LAST such cutoff avoids
    # assuming cutoff-specific thresholds are monotonic.
    upper_candidates = np.flatnonzero(
        cumulative <= upper_thresholds
    )

    upper = (
        int(upper_candidates[-1]) + 2
        if upper_candidates.size
        else 1
    )

    assert 1 <= lower <= upper <= D_MAX

    return lower, upper


rng = np.random.default_rng(20261102)
records = []

for exposure_name, exposure in exposures.items():
    assert np.all(exposure > 0)

    for scenario, intensity in intensities.items():
        truth = truths[scenario]
        competitor = 28 if scenario == "narrowing" else 9

        for depth in (1.0, 0.01, 0.001):
            expected = exposure * intensity * depth

            for replicate in range(N_REPLICATES):
                counts = rng.poisson(expected)

                global_interval = global_bounds(
                    counts, exposure
                )

                local_interval = cutoff_specific_bounds(
                    counts, exposure
                )

                assert (
                    (global_interval is None)
                    == (local_interval is None)
                )

                record = {
                    "exposure": exposure_name,
                    "scenario": scenario,
                    "depth": depth,
                    "replicate": replicate,
                    "contacts": int(counts.sum()),
                    "valid": local_interval is not None,
                    "global_covered": False,
                    "local_covered": False,
                    "global_discriminates": False,
                    "local_discriminates": False,
                    "global_span": np.nan,
                    "local_span": np.nan,
                }

                if local_interval is not None:
                    gl, gu = global_interval
                    ll, lu = local_interval

                    global_covered = gl <= truth <= gu
                    local_covered = ll <= truth <= lu

                    record.update({
                        "global_covered": global_covered,
                        "local_covered": local_covered,
                        "global_discriminates": (
                            global_covered
                            and not (gl <= competitor <= gu)
                        ),
                        "local_discriminates": (
                            local_covered
                            and not (ll <= competitor <= lu)
                        ),
                        "global_span": gu - gl,
                        "local_span": lu - ll,
                    })

                records.append(record)

df = pd.DataFrame(records)

assert len(df) == (
    len(exposures)
    * len(intensities)
    * 3
    * N_REPLICATES
)

output = Path(
    "research_artifacts/stage1/"
    "cutoff_exposure_rank_comparison.csv"
)

output.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(output, index=False)

print("GLOBAL VS CUTOFF-SPECIFIC EXPOSURE BOUNDS")
print("=" * 115)

print(
    f"{'Exposure':<13} {'Scenario':<15} {'Depth':>8} "
    f"{'Valid':>9} {'Global cov':>11} {'Local cov':>10} "
    f"{'Global disc':>12} {'Local disc':>11} "
    f"{'Global span':>12} {'Local span':>11}"
)
print("-" * 115)

for (exposure, scenario, depth), group in df.groupby(
    ["exposure", "scenario", "depth"],
    sort=False,
):
    valid = group[group.valid]

    def coverage(column):
        return (
            f"{100 * valid[column].mean():.1f}%"
            if len(valid) else "NA"
        )

    def span(column):
        return (
            f"{valid[column].median():.2f}"
            if len(valid) else "NA"
        )

    print(
        f"{exposure:<13} "
        f"{scenario:<15} "
        f"{depth:8.4f} "
        f"{len(valid):3d}/{len(group):<5d} "
        f"{coverage('global_covered'):>11} "
        f"{coverage('local_covered'):>10} "
        f"{100*group.global_discriminates.mean():11.1f}% "
        f"{100*group.local_discriminates.mean():10.1f}% "
        f"{span('global_span'):>12} "
        f"{span('local_span'):>11}"
    )

print()
print(f"Output: {output}")
print("Coverage is conditional on a valid interval.")
print("Discrimination is measured across all attempts.")
print("Intervals are deterministic functions of observed counts.")
print("No biological validation or two-window significance test.")
print("PASS: Cutoff-specific exposure comparison completed")
