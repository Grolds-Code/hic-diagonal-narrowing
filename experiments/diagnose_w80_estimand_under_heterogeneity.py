"""Distinguish population W80 from sample-specific latent W80.

Controlled simulations:
- Independent Poisson counts.
- Fixed heterogeneous genomic-pair intensities.
- Redrawn Gamma-Poisson pair effects.
- Gamma effects shared across varying numbers of distance bands.

A conditional Poisson interval targets the realized latent
profile when random effects are present. It does not
automatically cover the unconditional population profile.

These simulations do not establish biological coverage.
"""

import json
from math import sqrt
from pathlib import Path

import numpy as np

from hic_diagonal_narrowing.poisson_uncertainty import (
    cutoff_exposure_w80_interval_from_bands,
)

SEED = 20261028
TRIALS = 600
K = 64
DMAX = 64

PROFILES = {
    "Region_A_like": (8.4, 2478),
    "Region_B_like": (7.5, 4695),
}

MODELS = [
    "independent_poisson",
    "fixed_pair_heterogeneity",
    "redrawn_pair_gamma",
    "block_1",
    "block_2",
    "block_4",
    "block_8",
    "block_16",
    "block_64",
]


def w80_from_intensity(values):
    values = np.asarray(values, dtype=float)

    if (
        values.shape != (DMAX,)
        or not np.all(np.isfinite(values))
        or np.any(values < 0)
        or values.sum() <= 0
    ):
        raise ValueError("Invalid population intensity")

    return int(
        np.searchsorted(
            np.cumsum(values),
            0.8 * values.sum(),
        ) + 1
    )


def make_mean_profile(tau, total):
    d = np.arange(1, DMAX + 1)
    shape = np.exp(-d / tau)

    mu = shape * total / (K * shape.sum())

    return mu


def wilson(successes, n):
    z = 1.959963984540054
    p = successes / n
    denom = 1 + z * z / n

    center = (p + z*z/(2*n)) / denom

    margin = (
        z / denom
        * sqrt(
            p*(1-p)/n + z*z/(4*n*n)
        )
    )

    return [center - margin, center + margin]


def latent_band_means(rng, mu, model, fixed_means):
    """Return conditional Poisson means for all distance bands."""

    baseline = K * mu

    if model == "independent_poisson":
        return baseline.copy()

    if model == "fixed_pair_heterogeneity":
        # Pair-specific means drawn once per simulation
        # scenario, then held fixed across datasets.
        return fixed_means.copy()

    if model == "redrawn_pair_gamma":
        # Sum of K independent Gamma(kappa, mu/kappa)
        # pair-level intensity effects.
        #
        # Redrawn in every simulated dataset.
        kappa = 0.5

        return rng.gamma(
            shape=K * kappa,
            scale=mu / kappa,
        )

    if model.startswith("block_"):
        block = int(model.split("_")[1])

        if DMAX % block:
            raise RuntimeError("Block must divide DMAX")

        # Unit-mean Gamma effects with variance 1/8.
        effects = rng.gamma(
            shape=8.0,
            scale=1 / 8.0,
            size=DMAX // block,
        )

        return baseline * np.repeat(effects, block)

    raise ValueError(f"Unknown model: {model}")


def run_case(mu, model, profile_index, model_index):
    rng = np.random.default_rng(
        np.random.SeedSequence([
            SEED,
            profile_index,
            model_index,
        ])
    )

    baseline = K * mu
    population_target = w80_from_intensity(baseline)

    # Large fixed heterogeneity across pairs, without
    # creating extra sampling variance across datasets.
    pair_means = rng.gamma(
        shape=0.5,
        scale=mu[None, :] / 0.5,
        size=(K, DMAX),
    )

    fixed_means = pair_means.sum(axis=0)

    exposures = np.full(DMAX, K, dtype=np.int64)

    covered_population = 0
    covered_latent = 0
    latent_shifted = 0

    population_misses_below = 0
    population_misses_above = 0

    interval_spans = []
    latent_distances = []
    realized_targets = []

    for _ in range(TRIALS):
        band_means = latent_band_means(
            rng,
            mu,
            model,
            fixed_means,
        )

        if (
            np.any(~np.isfinite(band_means))
            or np.any(band_means < 0)
        ):
            raise RuntimeError("Invalid latent intensity")

        latent_target = w80_from_intensity(
            band_means
        )

        counts = rng.poisson(band_means).astype(
            np.int64
        )

        if counts.sum() == 0:
            raise RuntimeError(
                "Unexpected zero-contact realization"
            )

        result = (
            cutoff_exposure_w80_interval_from_bands(
                counts,
                exposures,
                confidence_level=0.95,
            )
        )

        lo = int(result["cutoff_lower"])
        hi = int(result["cutoff_upper"])
        observed = int(result["observed_w80"])

        if not 1 <= lo <= observed <= hi <= DMAX:
            raise RuntimeError("Invalid W80 interval")

        covered_population += int(
            lo <= population_target <= hi
        )

        covered_latent += int(
            lo <= latent_target <= hi
        )

        latent_shifted += int(
            latent_target != population_target
        )

        population_misses_below += int(
            population_target < lo
        )

        population_misses_above += int(
            population_target > hi
        )

        interval_spans.append(hi - lo)
        latent_distances.append(
            abs(latent_target - population_target)
        )
        realized_targets.append(latent_target)

    if (
        covered_population
        + population_misses_below
        + population_misses_above
        != TRIALS
    ):
        raise RuntimeError("Coverage accounting mismatch")

    return {
        "model": model,
        "trials": TRIALS,
        "population_target_w80": population_target,
        "coverage_population_target": (
            covered_population / TRIALS
        ),
        "population_coverage_wilson_95": wilson(
            covered_population, TRIALS
        ),
        "coverage_latent_target": (
            covered_latent / TRIALS
        ),
        "latent_coverage_wilson_95": wilson(
            covered_latent, TRIALS
        ),
        "fraction_latent_target_shifted": (
            latent_shifted / TRIALS
        ),
        "median_absolute_latent_shift_bins": float(
            np.median(latent_distances)
        ),
        "median_interval_span_bins": float(
            np.median(interval_spans)
        ),
        "population_target_miss_below": (
            population_misses_below / TRIALS
        ),
        "population_target_miss_above": (
            population_misses_above / TRIALS
        ),
        "range_latent_w80": [
            int(min(realized_targets)),
            int(max(realized_targets)),
        ],
        "fixed_heterogeneity_target_w80": (
            w80_from_intensity(fixed_means)
            if model == "fixed_pair_heterogeneity"
            else None
        ),
    }


def main():
    output = Path(
        "research_artifacts/stage1/"
        "w80_estimand_heterogeneity_diagnostic_v1.json"
    )

    if output.exists():
        raise RuntimeError(
            "Output already exists; refusing to overwrite"
        )

    results = {}

    print("\nW80 ESTIMAND AND HETEROGENEITY DIAGNOSTIC")
    print("=" * 79)
    print("Nominal interval confidence: 95%")
    print("Trials per model and profile:", TRIALS)

    for pi, (name, (tau, total)) in enumerate(
        PROFILES.items()
    ):
        mu = make_mean_profile(tau, total)

        target = w80_from_intensity(mu)

        results[name] = {
            "unconditional_mean_w80": target,
            "expected_total_contacts": total,
            "models": {},
        }

        print(f"\n{name}: population W80={target} bins")
        print(
            f"{'Model':29s} "
            f"{'Population':>11s} "
            f"{'Latent':>10s} "
            f"{'Target shifted':>15s} "
            f"{'Median span':>12s}"
        )

        for mi, model in enumerate(MODELS):
            result = run_case(
                mu, model, pi, mi
            )

            results[name]["models"][model] = result

            print(
                f"{model:29s} "
                f"{result['coverage_population_target']:>10.1%} "
                f"{result['coverage_latent_target']:>10.1%} "
                f"{result['fraction_latent_target_shifted']:>14.1%} "
                f"{result['median_interval_span_bins']:>11.1f}"
            )

    payload = {
        "experiment": "w80_estimand_under_heterogeneity",
        "seed": SEED,
        "trials_per_scenario": TRIALS,
        "confidence_level": 0.95,
        "exposures_per_distance": K,
        "distance_bands": DMAX,
        "fixed_pair_effect_shape": 0.5,
        "redrawn_pair_effect_shape": 0.5,
        "block_effect_shape": 8.0,
        "results": results,
        "interpretation": (
            "Coverage of unconditional population W80 versus "
            "realized sample-specific latent W80. Synthetic "
            "random-effect models are hypothetical and not "
            "fitted to HCT116 biological replicates."
        ),
    }

    output.parent.mkdir(parents=True, exist_ok=True)

    output.write_text(
        json.dumps(payload, indent=2) + "\n"
    )

    print("\nSaved:", output)
    print(
        "These are synthetic repeated-sampling diagnostics, "
        "not HCT116 biological coverage estimates."
    )


if __name__ == "__main__":
    main()
