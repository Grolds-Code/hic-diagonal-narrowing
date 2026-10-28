"""Audit randomization sensitivity of equal-exposure W80 intervals.

Each synthetic observation is held fixed while the equal-exposure
thinning is repeated independently.

This measures algorithmic repeatability, not population coverage.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binom

reference = pd.read_csv(
    "research_artifacts/stage1/"
    "poisson_unequal_support_calibration.csv"
).set_index(["scenario", "depth", "replicate"])

N_BINS = 256
D_MAX = 64
N_REPEATS = 100
SELECTED_REPLICATES = set(range(0, 200, 10))

SCENARIOS = ("baseline", "narrowing", "coverage_drop")
DEPTHS = (1.0, 0.01, 0.001, 0.0001)

distances = np.arange(1, D_MAX + 1)

mask_rng = np.random.default_rng(20261028)
valid_bins = np.ones(N_BINS, dtype=bool)
valid_bins[
    mask_rng.choice(N_BINS, size=58, replace=False)
] = False

eligible = np.zeros(D_MAX, dtype=np.int64)

for d in distances:
    rows = np.arange(N_BINS - d)
    midpoint = rows + d / 2

    selected = (
        (midpoint >= 96)
        & (midpoint < 160)
        & valid_bins[rows]
        & valid_bins[rows + d]
    )

    eligible[d - 1] = selected.sum()

assert eligible.min() == 36
assert eligible.max() == 44

retention = eligible.min() / eligible

baseline = 90 * np.exp(-distances / 18)
narrowing = 90 * np.exp(-distances / 5)

narrowing *= (
    np.sum(eligible * baseline)
    / np.sum(eligible * narrowing)
)

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


def rank_bounds(counts):
    n = int(counts.sum())

    if n == 0:
        return None

    lo_rank = int(binom.ppf(0.025, n, 0.80))
    hi_rank = int(binom.isf(0.025, n, 0.80)) + 1

    cumulative = np.cumsum(counts)

    lo = (
        1 if lo_rank < 1
        else int(np.searchsorted(cumulative, lo_rank) + 1)
    )

    hi = (
        D_MAX if hi_rank > n
        else int(np.searchsorted(cumulative, hi_rank) + 1)
    )

    assert 1 <= lo <= hi <= D_MAX
    return lo, hi


count_rng = np.random.default_rng(810241)

case_records = []
repeat_records = []

for si, scenario in enumerate(SCENARIOS):
    truth = truths[scenario]
    alternative = 28 if scenario == "narrowing" else 9

    for di, depth in enumerate(DEPTHS):
        expected = eligible * intensities[scenario] * depth

        for rep in range(200):
            counts = count_rng.poisson(expected).astype(np.int64)

            expected_contacts = int(
                reference.loc[
                    (scenario, depth, rep), "contacts"
                ]
            )

            assert int(counts.sum()) == expected_contacts

            if rep not in SELECTED_REPLICATES:
                continue

            intervals = set()
            decisions = set()
            valid_count = 0
            retained_counts = []
            spans = []

            for run in range(N_REPEATS):
                rng = np.random.default_rng(
                    20261030
                    + si * 10000000
                    + di * 100000
                    + rep * 1000
                    + run
                )

                thinned = rng.binomial(
                    counts, retention
                )

                retained = int(thinned.sum())
                retained_counts.append(retained)

                bounds = rank_bounds(thinned)

                if bounds is None:
                    state = "abstain"
                    lower = np.nan
                    upper = np.nan
                else:
                    valid_count += 1
                    lower, upper = bounds

                    intervals.add(bounds)
                    spans.append(upper - lower)

                    discriminates = (
                        lower <= truth <= upper
                        and not (
                            lower <= alternative <= upper
                        )
                    )

                    state = (
                        "discriminate"
                        if discriminates
                        else "not_discriminate"
                    )

                decisions.add(state)

                repeat_records.append({
                    "scenario": scenario,
                    "depth": depth,
                    "replicate": rep,
                    "thinning_run": run,
                    "original_contacts": int(counts.sum()),
                    "retained_contacts": retained,
                    "lower": lower,
                    "upper": upper,
                    "decision": state,
                })

            valid_decisions = decisions - {"abstain"}

            case_records.append({
                "scenario": scenario,
                "depth": depth,
                "replicate": rep,
                "valid_fraction": valid_count / N_REPEATS,
                "unique_intervals": len(intervals),
                "decision_switch": len(valid_decisions) > 1,
                "abstention_switch": (
                    "abstain" in decisions
                    and valid_count > 0
                ),
                "median_retained": float(
                    np.median(retained_counts)
                ),
                "median_span": (
                    float(np.median(spans))
                    if spans else np.nan
                ),
            })

cases = pd.DataFrame(case_records)
repeats = pd.DataFrame(repeat_records)

assert len(cases) == 240
assert len(repeats) == 24000

output = Path(
    "research_artifacts/stage1/"
    "equalized_rank_repeatability.csv"
)

output.parent.mkdir(parents=True, exist_ok=True)
repeats.to_csv(output, index=False)

print("EQUAL-EXPOSURE RANK INTERVAL REPEATABILITY")
print("=" * 103)
print(f"Fixed synthetic observations: {len(cases)}")
print(f"Thinning repeats per observation: {N_REPEATS}")
print(f"Total interval attempts: {len(repeats)}")
print()

print(
    f"{'Scenario':<15} {'Depth':>8} "
    f"{'Cases':>6} {'Valid %':>9} "
    f"{'Med unique':>11} {'Decision flips':>15} "
    f"{'Abstain flips':>14} {'Med retained':>13}"
)
print("-" * 103)

for (scenario, depth), group in cases.groupby(
    ["scenario", "depth"], sort=False
):
    print(
        f"{scenario:<15} {depth:8.4f} "
        f"{len(group):6d} "
        f"{100*group.valid_fraction.mean():8.1f}% "
        f"{group.unique_intervals.median():11.1f} "
        f"{int(group.decision_switch.sum()):6d}/{len(group):<8d} "
        f"{int(group.abstention_switch.sum()):6d}/{len(group):<7d} "
        f"{group.median_retained.median():13.1f}"
    )

print()
print(f"Output: {output}")
print("Decision flips compare repeated thinning of the SAME data.")
print("This is not an independent population coverage experiment.")
print("PASS: Equal-exposure interval repeatability audited")
