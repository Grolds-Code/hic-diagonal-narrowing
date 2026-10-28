"""Stress-test W80 confidence intervals under model misspecification.

All simulations have the same deterministic mean contact profile
within each synthetic case.

The target W80 is the W80 of the unconditional mean profile.

Results measure synthetic repeated-sampling coverage, NOT
biological confidence coverage for real Hi-C data.
"""

import json
from math import sqrt
from pathlib import Path

import numpy as np

from hic_diagonal_narrowing.poisson_uncertainty import (
    cutoff_exposure_w80_interval_from_bands,
)

SEED = 20261028
TRIALS = 1000
DMAX = 64
EXPOSURE = 64
CONFIDENCE = 0.95

# Synthetic exponential profiles approximating the observed
# full-depth W80 and contact totals of the two candidate windows.
PROFILES = {
    "Region_A_like": {
        "tau": 8.4,
        "total_expected_contacts": 2478,
        "expected_w80": 14,
    },
    "Region_B_like": {
        "tau": 7.5,
        "total_expected_contacts": 4695,
        "expected_w80": 13,
    },
}

MODELS = (
    "independent_poisson",
    "global_depth_gamma",
    "independent_pair_gamma_poisson",
    "correlated_distance_block_gamma",
)


def wilson_interval(successes, n):
    """95% Wilson interval for a simulated coverage proportion."""
    z = 1.959963984540054
    p = successes / n
    denominator = 1 + z * z / n
    center = (p + z * z / (2 * n)) / denominator
    half_width = (
        z / denominator
        * sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    )
    return [center - half_width, center + half_width]


def make_profile(tau, total):
    distances = np.arange(1, DMAX + 1)
    shape = np.exp(-distances / tau)

    # mu[d] is expected contact count per eligible pair.
    mu = shape * total / (EXPOSURE * shape.sum())

    assert np.isclose(EXPOSURE * mu.sum(), total)

    return mu


def generate_band_counts(rng, mu, model):
    band_mean = EXPOSURE * mu

    if model == "independent_poisson":
        # Y_d ~ Poisson(K * mu_d), independent across d.
        return rng.poisson(band_mean)

    if model == "global_depth_gamma":
        # One random multiplier shared by all bands.
        # Marginal counts are correlated, but conditional on
        # N=sum(Y_d), the probabilities remain proportional
        # to K * mu_d.
        global_multiplier = rng.gamma(
            shape=4.0, scale=1 / 4.0
        )
        return rng.poisson(
            band_mean * global_multiplier
        )

    if model == "independent_pair_gamma_poisson":
        # Each of K contact pairs has an independent
        # Gamma-Poisson random effect, redrawn per dataset.
        #
        # The aggregate is negative binomial with:
        # mean = K * mu_d
        # variance = K * mu_d + K * mu_d^2 / kappa
        kappa = 0.5

        return rng.negative_binomial(
            n=EXPOSURE * kappa,
            p=kappa / (kappa + mu),
        )

    if model == "correlated_distance_block_gamma":
        # Eight blocks of eight adjacent distances.
        # Each block shares a random intensity multiplier.
        # Multipliers are redrawn in each synthetic dataset.
        block_effects = rng.gamma(
            shape=8.0,
            scale=1 / 8.0,
            size=DMAX // 8,
        )

        multipliers = np.repeat(block_effects, 8)

        return rng.poisson(
            band_mean * multipliers
        )

    raise ValueError(f"Unknown model: {model}")


def simulate_case(rng, mu, true_w80, model):
    exposures = np.full(
        DMAX, EXPOSURE, dtype=np.int64
    )

    covered = 0
    lower_misses = 0
    upper_misses = 0
    empty = 0

    widths = []
    errors = []
    correct = 0

    for _ in range(TRIALS):
        counts = np.asarray(
            generate_band_counts(rng, mu, model),
            dtype=np.int64,
        )

        if counts.shape != (DMAX,) or np.any(counts < 0):
            raise RuntimeError("Invalid simulated counts")

        if counts.sum() == 0:
            empty += 1
            continue

        result = (
            cutoff_exposure_w80_interval_from_bands(
                counts,
                exposures,
                confidence_level=CONFIDENCE,
            )
        )

        observed = int(result["observed_w80"])
        lower = int(result["cutoff_lower"])
        upper = int(result["cutoff_upper"])

        if not 1 <= lower <= observed <= upper <= DMAX:
            raise RuntimeError("Invalid returned interval")

        widths.append(upper - lower)
        errors.append(abs(observed - true_w80))
        correct += int(observed == true_w80)

        covered += int(lower <= true_w80 <= upper)
        lower_misses += int(true_w80 < lower)
        upper_misses += int(true_w80 > upper)

    n = TRIALS - empty

    if n < TRIALS * 0.99:
        raise RuntimeError(
            "Too many empty datasets for reliable comparison"
        )

    if covered + lower_misses + upper_misses != n:
        raise RuntimeError("Coverage accounting mismatch")

    return {
        "model": model,
        "trials_requested": TRIALS,
        "valid_trials": n,
        "empty_trials": empty,
        "true_population_w80_bins": true_w80,
        "empirical_coverage": covered / n,
        "coverage_wilson_95": wilson_interval(
            covered, n
        ),
        "lower_tail_miss_fraction": lower_misses / n,
        "upper_tail_miss_fraction": upper_misses / n,
        "median_interval_span_bins": float(
            np.median(widths)
        ),
        "median_absolute_w80_error_bins": float(
            np.median(errors)
        ),
        "fraction_point_estimate_exact": correct / n,
    }


def main():
    rng = np.random.default_rng(SEED)
    results = {}

    print("\nW80 MODEL-MISSPECIFICATION STRESS TEST")
    print("=" * 72)
    print("Nominal confidence: 95%")
    print("Simulated datasets per case:", TRIALS)
    print("Eligible pairs per distance:", EXPOSURE)

    for name, config in PROFILES.items():
        mu = make_profile(
            config["tau"],
            config["total_expected_contacts"],
        )

        true_w80 = int(
            np.searchsorted(
                np.cumsum(mu),
                0.80 * mu.sum(),
            ) + 1
        )

        if true_w80 != config["expected_w80"]:
            raise RuntimeError(
                f"Unexpected population W80 for {name}"
            )

        print(f"\n{name}")
        print(
            f"Population W80: {true_w80} bins; "
            f"expected contacts: "
            f"{config['total_expected_contacts']}"
        )

        print(
            f"{'Model':37s} {'Coverage':>9s} "
            f"{'Wilson 95%':>17s} {'Median span':>12s}"
        )

        results[name] = {
            "profile": config,
            "true_population_w80_bins": true_w80,
            "models": {},
        }

        for model in MODELS:
            result = simulate_case(
                rng, mu, true_w80, model
            )

            results[name]["models"][model] = result

            lo, hi = result["coverage_wilson_95"]

            print(
                f"{model:37s} "
                f"{result['empirical_coverage']:>8.1%} "
                f"[{lo:.1%}, {hi:.1%}] "
                f"{result['median_interval_span_bins']:>11.1f}"
            )

    output = Path(
        "research_artifacts/stage1/"
        "w80_model_misspecification_stress_v1.json"
    )

    output.parent.mkdir(parents=True, exist_ok=True)

    if output.exists():
        raise RuntimeError(
            "Output already exists; refusing to overwrite"
        )

    payload = {
        "experiment": "w80_model_misspecification_stress",
        "seed": SEED,
        "trials_per_case": TRIALS,
        "nominal_confidence": CONFIDENCE,
        "exposure_per_band": EXPOSURE,
        "distance_bands": DMAX,
        "models": {
            "independent_poisson": (
                "Independent Poisson aggregate band counts"
            ),
            "global_depth_gamma": (
                "Shared Gamma(shape=4) depth multiplier "
                "preserving conditional multinomial probabilities"
            ),
            "independent_pair_gamma_poisson": (
                "Independent pair-level Gamma-Poisson effects; "
                "negative-binomial shape kappa=0.5"
            ),
            "correlated_distance_block_gamma": (
                "Shared Gamma(shape=8) multiplier within "
                "each eight-distance block"
            ),
        },
        "results": results,
        "interpretation": (
            "Synthetic coverage for W80 of unconditional "
            "population mean profiles. These models are "
            "hypothetical stress scenarios, not empirically "
            "fitted HCT116 count-generation models."
        ),
    }

    output.write_text(
        json.dumps(payload, indent=2) + "\n"
    )

    print("\nSaved:", output)
    print(
        "No inference about biological confidence coverage "
        "is made by this experiment."
    )


if __name__ == "__main__":
    main()
