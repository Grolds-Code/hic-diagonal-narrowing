"""Verify reusable cutoff-exposure W80 intervals against the
previous experimental implementation and saved calibration.

Reproduces all 7,200 synthetic observations using original seeds.
"""

import ast
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binom

from hic_diagonal_narrowing.poisson_uncertainty import (
    cutoff_exposure_w80_interval_from_bands,
)

reference = pd.read_csv(
    "research_artifacts/stage1/"
    "cutoff_exposure_rank_comparison.csv"
)

assert len(reference) == 7200

keys = ["exposure", "scenario", "depth", "replicate"]
assert not reference.duplicated(keys).any()
reference = reference.set_index(keys)

# Load only the original experimental function,
# without executing its simulation.
script = Path(
    "experiments/compare_cutoff_exposure_rank.py"
)

tree = ast.parse(script.read_text())

functions = [
    node for node in tree.body
    if isinstance(node, ast.FunctionDef)
    and node.name == "cutoff_specific_bounds"
]

assert len(functions) == 1

namespace = {
    "np": np,
    "binom": binom,
    "D_MAX": 64,
    "ALPHA": 0.05,
    "TARGET": 0.80,
}

exec(
    compile(
        ast.Module(body=functions, type_ignores=[]),
        str(script),
        "exec",
    ),
    namespace,
)

original_function = namespace["cutoff_specific_bounds"]

distances = np.arange(1, 65)

mask_rng = np.random.default_rng(20261028)
valid_bins = np.ones(256, dtype=bool)
valid_bins[
    mask_rng.choice(256, size=58, replace=False)
] = False

masked_exposure = np.zeros(64, dtype=np.int64)

for d in distances:
    rows = np.arange(256 - d)
    midpoint = rows + d / 2

    selected = (
        (midpoint >= 96)
        & (midpoint < 160)
        & valid_bins[rows]
        & valid_bins[rows + d]
    )

    masked_exposure[d - 1] = selected.sum()

exposures = {
    "equal": np.full(64, 40, dtype=np.int64),
    "masked": masked_exposure,
    "decreasing": np.rint(
        np.linspace(70, 10, 64)
    ).astype(np.int64),
    "increasing": np.rint(
        np.linspace(10, 70, 64)
    ).astype(np.int64),
}

intensities = {
    "baseline": 90 * np.exp(-distances / 18),
    "narrowing": 90 * np.exp(-distances / 5),
    "coverage_drop": 22.5 * np.exp(-distances / 18),
}

truths = {
    "baseline": 28,
    "narrowing": 9,
    "coverage_drop": 28,
}

rng = np.random.default_rng(20261102)

checked = 0
valid = 0
abstentions = 0
mismatches = []

for exposure_name, exposure in exposures.items():
    for scenario, intensity in intensities.items():
        truth = truths[scenario]
        competitor = 28 if scenario == "narrowing" else 9

        for depth in (1.0, 0.01, 0.001):
            expected_counts = exposure * intensity * depth

            for replicate in range(200):
                counts = rng.poisson(expected_counts)

                key = (
                    exposure_name,
                    scenario,
                    depth,
                    replicate,
                )

                record = reference.loc[key]

                assert int(counts.sum()) == int(
                    record["contacts"]
                )

                original = original_function(
                    counts, exposure
                )

                if original is None:
                    assert not bool(record["valid"])
                    abstentions += 1

                    try:
                        cutoff_exposure_w80_interval_from_bands(
                            counts, exposure
                        )
                    except ValueError as exc:
                        assert "No informative" in str(exc)
                    else:
                        raise AssertionError(
                            f"Zero-contact trial accepted: {key}"
                        )

                else:
                    assert bool(record["valid"])

                    result = (
                        cutoff_exposure_w80_interval_from_bands(
                            counts,
                            exposure,
                            confidence_level=0.95,
                        )
                    )

                    actual = (
                        result["cutoff_lower"],
                        result["cutoff_upper"],
                    )

                    lower, upper = actual

                    covered = lower <= truth <= upper

                    discriminates = (
                        covered
                        and not (
                            lower <= competitor <= upper
                        )
                    )

                    matches = (
                        actual == original
                        and covered == bool(
                            record["local_covered"]
                        )
                        and discriminates == bool(
                            record["local_discriminates"]
                        )
                        and (upper - lower) == int(
                            record["local_span"]
                        )
                    )

                    if not matches:
                        mismatches.append({
                            "trial": key,
                            "original": original,
                            "actual": actual,
                        })

                    valid += 1

                checked += 1

assert checked == 7200
assert valid + abstentions == checked

print("CUTOFF-EXPOSURE IMPLEMENTATION VERIFICATION")
print("=" * 72)
print(f"Total synthetic trials checked: {checked}")
print(f"Valid intervals checked: {valid}")
print(f"Zero-contact abstentions verified: {abstentions}")
print(f"Mismatches: {len(mismatches)}")

if mismatches:
    for mismatch in mismatches[:10]:
        print(mismatch)

    raise AssertionError(
        "Reusable implementation differs from calibration"
    )

print("PASS: Interval endpoints match previous implementation")
print("PASS: Coverage and discrimination outcomes match")
print("PASS: Interval spans match saved calibration")
print("PASS: All 7,200 trials verified")
