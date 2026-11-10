"""Controlled stress test of exploratory spatial-coverage QC rules."""

from pathlib import Path

import numpy as np
import pandas as pd

from hic_diagonal_narrowing.synthetic import simulate_contact_map
from hic_diagonal_narrowing.width import estimate_w80

N = 256
REGION = (96, 160)
DMAX = 64
REPLICATES = 40
DEPTHS = (0.03, 0.01)
SCORE_THRESHOLD = 0.20

ROWS, COLS = np.triu_indices(N, k=1)
NEAR = (COLS - ROWS) <= DMAX

OUTPUT = Path(
    "research_artifacts/stage1/"
    "stage1_qc_dropout_stress_v1.csv"
)

if OUTPUT.exists():
    raise RuntimeError(f"Refusing to overwrite: {OUTPUT}")


def thin(matrix, depth, seed):
    rng = np.random.default_rng(seed)
    sampled = rng.binomial(
        matrix[ROWS, COLS].astype(np.int64),
        depth,
    )

    result = np.zeros_like(matrix)
    result[ROWS, COLS] = sampled
    result[COLS, ROWS] = sampled
    return result


def inject_gap(matrix, start, length, retention, seed):
    """Reduce counts for pairs touching specified genomic bins."""
    result = matrix.copy()

    affected = (
        ((ROWS >= start) & (ROWS < start + length))
        | ((COLS >= start) & (COLS < start + length))
    )

    if retention == 0:
        values = np.zeros(
            int(affected.sum()), dtype=np.int64
        )
    else:
        rng = np.random.default_rng(seed)
        values = rng.binomial(
            result[ROWS[affected], COLS[affected]]
            .astype(np.int64),
            retention,
        )

    result[ROWS[affected], COLS[affected]] = values
    result[COLS[affected], ROWS[affected]] = values

    return result


def longest_run(mask):
    padded = np.r_[False, mask, False].astype(int)
    transitions = np.diff(padded)

    starts = np.flatnonzero(transitions == 1)
    ends = np.flatnonzero(transitions == -1)

    return (
        int(np.max(ends - starts))
        if len(starts) else 0
    )


def assess(matrix, reference_width):
    width = estimate_w80(
        matrix, REGION,
        max_separation_bins=DMAX,
    )

    score = float(
        np.log2(reference_width / width)
    )

    observed = (
        NEAR & (matrix[ROWS, COLS] > 0)
    )

    partners = np.bincount(
        np.r_[
            ROWS[observed],
            COLS[observed],
        ],
        minlength=N,
    )

    # The 64-kb endpoint span surrounding the
    # central 64-bin measurement window.
    endpoint = partners[64:192]

    if len(endpoint) != 128:
        raise RuntimeError("Unexpected endpoint span")

    median = float(np.median(endpoint))

    if median <= 0:
        raise RuntimeError("No endpoint contact support")

    zero_run = longest_run(endpoint == 0)
    low10_run = longest_run(
        endpoint < 0.10 * median
    )
    low25_run = longest_run(
        endpoint < 0.25 * median
    )

    warning_a = zero_run >= 2
    warning_b = low10_run >= 2
    warning_c = low25_run >= 2

    raw_call = score > SCORE_THRESHOLD

    return {
        "w80_bins": width,
        "narrowing_score": score,
        "raw_call": raw_call,
        "warning_A": warning_a,
        "warning_B": warning_b,
        "warning_C": warning_c,
        "accepted_call_A": raw_call and not warning_a,
        "accepted_call_B": raw_call and not warning_b,
        "accepted_call_C": raw_call and not warning_c,
        "longest_zero_run": zero_run,
        "longest_low10_run": low10_run,
        "longest_low25_run": low25_run,
    }


records = []

for replicate in range(REPLICATES):
    reference_full = simulate_contact_map(
        "baseline", seed=900000 + replicate
    )

    baseline_full = simulate_contact_map(
        "baseline", seed=100000 + replicate
    )

    narrowing_full = simulate_contact_map(
        "narrowing", seed=200000 + replicate
    )

    coverage_full = simulate_contact_map(
        "coverage_drop", seed=300000 + replicate
    )

    # Vary the location independently of the QC outcomes.
    gap_start = (83, 102, 128, 151, 173)[replicate % 5]

    for depth_index, depth in enumerate(DEPTHS):
        prefix = replicate * 100 + depth_index * 10

        reference = thin(
            reference_full, depth, 5000000 + prefix
        )
        baseline = thin(
            baseline_full, depth, 6000000 + prefix
        )
        narrowing = thin(
            narrowing_full, depth, 7000000 + prefix
        )
        uniform = thin(
            coverage_full, depth, 8000000 + prefix
        )

        reference_width = estimate_w80(
            reference, REGION,
            max_separation_bins=DMAX,
        )

        cases = {
            "baseline": baseline,
            "true_narrowing": narrowing,
            "uniform_coverage_loss": uniform,
            "single_zero_bin": inject_gap(
                baseline, gap_start, 1, 0,
                9000000 + prefix,
            ),
            "four_zero_bins": inject_gap(
                baseline, gap_start, 4, 0,
                9100000 + prefix,
            ),
            "four_90pct_loss": inject_gap(
                baseline, gap_start, 4, 0.10,
                9200000 + prefix,
            ),
            "narrowing_plus_four_zeros": inject_gap(
                narrowing, gap_start, 4, 0,
                9300000 + prefix,
            ),
        }

        for scenario, matrix in cases.items():
            metrics = assess(
                matrix, reference_width
            )

            records.append({
                "replicate": replicate,
                "depth": depth,
                "scenario": scenario,
                "gap_start_bin": (
                    gap_start if scenario not in (
                        "baseline",
                        "true_narrowing",
                        "uniform_coverage_loss",
                    ) else -1
                ),
                "reference_w80_bins": reference_width,
                **metrics,
            })

result = pd.DataFrame(records)

expected = REPLICATES * len(DEPTHS) * 7

if len(result) != expected:
    raise RuntimeError("Incomplete stress test")

if result.isna().any().any():
    raise RuntimeError("Unexpected missing measurements")

summary = (
    result.groupby(
        ["depth", "scenario"], sort=False
    )
    .agg(
        trials=("raw_call", "size"),
        raw_call_rate=("raw_call", "mean"),
        flag_A=("warning_A", "mean"),
        flag_B=("warning_B", "mean"),
        flag_C=("warning_C", "mean"),
        retained_A=("accepted_call_A", "mean"),
        retained_B=("accepted_call_B", "mean"),
        retained_C=("accepted_call_C", "mean"),
    )
    .reset_index()
)

if not summary["trials"].eq(REPLICATES).all():
    raise RuntimeError("Unequal scenario sample sizes")

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
result.to_csv(OUTPUT, index=False)

print("\nSTAGE 1 — SYNTHETIC COVERAGE-QC STRESS TEST")
print("=" * 100)
print("Simulation replicates:", REPLICATES)
print("Sequencing-depth fractions:", DEPTHS)
print("Scenarios:", result["scenario"].nunique())
print("Total simulation trials:", len(result))
print("Narrowing-score threshold: > 0.20")
print("Window geometry: 64 bins, 64 maximum separation")
print("QC endpoint span: 128 bins")
print("\nRates below are proportions from 0 to 1.\n")

print(summary.to_string(
    index=False,
    float_format=lambda x: f"{x:.3f}",
))

print("\nIMPORTANT LIMITATIONS")
print("- The synthetic model is not biological ground truth.")
print("- True-narrowing simulations use a strong engineered signal.")
print("- The reference is an independently simulated baseline.")
print("- This tests single-window behavior, not full-scan selection.")
print("- Genomic mappability and copy-number variation are not modeled.")
print("- QC thresholds are exploratory and not calibrated.")
print("- No existing Stage 1 calls have been altered.")
print("Saved:", OUTPUT)
print("\nStep 382: PASS")
