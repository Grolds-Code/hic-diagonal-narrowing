"""Cluster-level uncertainty for fixed heterogeneous W80 maps.

Independent sampling units are intensity maps (n=40),
not the repeated Poisson observations within each map.

Bootstrap intervals describe uncertainty across the
specified synthetic map-generating distribution.
"""

import json
from pathlib import Path

import numpy as np

ROOT = Path("research_artifacts/stage1")

SOURCE = ROOT / "w80_fixed_heterogeneity_map_audit_v1.json"
OUTPUT = ROOT / "w80_map_level_uncertainty_v1.json"

SEED = 20261029
BOOTSTRAPS = 10000


def bootstrap_coverage(maps, rng):
    population = np.asarray(
        [m["reference_population_inclusion"] for m in maps],
        dtype=float,
    )

    latent = np.asarray(
        [m["conditional_latent_inclusion"] for m in maps],
        dtype=float,
    )

    shifted = np.asarray(
        [m["latent_target_shift_bins"] != 0 for m in maps],
        dtype=bool,
    )

    n = len(maps)

    if n != 40:
        raise RuntimeError("Expected 40 independent maps")

    if (
        not np.all(np.isfinite(population))
        or not np.all(np.isfinite(latent))
        or np.any(population < 0)
        or np.any(population > 1)
        or np.any(latent < 0)
        or np.any(latent > 1)
    ):
        raise RuntimeError("Invalid per-map coverage")

    ids = [m["map_id"] for m in maps]

    if len(set(ids)) != n:
        raise RuntimeError("Duplicate map identifiers")

    # Resample entire maps, not individual Poisson trials.
    samples = rng.integers(
        low=0,
        high=n,
        size=(BOOTSTRAPS, n),
    )

    population_boot = population[samples].mean(axis=1)
    latent_boot = latent[samples].mean(axis=1)
    difference_boot = latent_boot - population_boot

    def interval(values):
        return [
            float(np.quantile(values, 0.025)),
            float(np.quantile(values, 0.975)),
        ]

    difference = latent - population

    return {
        "independent_map_count": n,
        "poisson_observations_per_map": 80,
        "coverage_population": float(population.mean()),
        "coverage_population_bootstrap_95": interval(
            population_boot
        ),
        "coverage_fixed_map": float(latent.mean()),
        "coverage_fixed_map_bootstrap_95": interval(
            latent_boot
        ),
        "coverage_difference_fixed_minus_population": float(
            difference.mean()
        ),
        "coverage_difference_bootstrap_95": interval(
            difference_boot
        ),
        "fraction_shifted_maps": float(shifted.mean()),
        "shifted_map_count": int(shifted.sum()),
        "population_coverage_shifted_maps": (
            float(population[shifted].mean())
            if shifted.any() else None
        ),
        "population_coverage_unshifted_maps": (
            float(population[~shifted].mean())
            if (~shifted).any() else None
        ),
        "fraction_bootstrap_population_coverage_below_95pct": (
            float(np.mean(population_boot < 0.95))
        ),
    }


def main():
    source = json.loads(SOURCE.read_text())

    if source["maps_per_profile"] != 40:
        raise RuntimeError("Unexpected map count")

    if source["replicates_per_map"] != 80:
        raise RuntimeError("Unexpected replicate count")

    if OUTPUT.exists():
        raise RuntimeError(
            "Output already exists; refusing to overwrite"
        )

    rng = np.random.default_rng(SEED)
    results = {}

    print("\nW80 MAP-LEVEL COVERAGE UNCERTAINTY")
    print("=" * 65)
    print("Independent maps per profile: 40")
    print("Poisson replicates per map: 80")
    print("Map-level bootstrap replicates:", BOOTSTRAPS)

    for name, profile in source["results"].items():
        maps = profile["maps"]

        result = bootstrap_coverage(maps, rng)

        if not np.isclose(
            result["coverage_population"],
            profile["mean_reference_population_inclusion"],
            atol=1e-12,
        ):
            raise RuntimeError("Population summary mismatch")

        if not np.isclose(
            result["coverage_fixed_map"],
            profile["mean_within_map_latent_inclusion"],
            atol=1e-12,
        ):
            raise RuntimeError("Latent summary mismatch")

        results[name] = result

        print(f"\n{name}")

        for label, key, ci_key in (
            (
                "Population-target coverage",
                "coverage_population",
                "coverage_population_bootstrap_95",
            ),
            (
                "Fixed-map-target coverage",
                "coverage_fixed_map",
                "coverage_fixed_map_bootstrap_95",
            ),
            (
                "Coverage difference",
                "coverage_difference_fixed_minus_population",
                "coverage_difference_bootstrap_95",
            ),
        ):
            value = result[key]
            lo, hi = result[ci_key]

            print(
                f"  {label}: "
                f"{value:.1%} "
                f"[{lo:.1%}, {hi:.1%}]"
            )

        print(
            "  Shifted maps:",
            result["shifted_map_count"],
            "/ 40",
        )

        print(
            "  Population coverage among shifted maps:",
            f"{result['population_coverage_shifted_maps']:.1%}"
        )

        print(
            "  Population coverage among unshifted maps:",
            f"{result['population_coverage_unshifted_maps']:.1%}"
        )

    OUTPUT.write_text(
        json.dumps({
            "experiment": "w80_map_level_uncertainty",
            "source_file": str(SOURCE),
            "bootstrap_seed": SEED,
            "bootstrap_replicates": BOOTSTRAPS,
            "independent_unit": "synthetic_intensity_map",
            "results": results,
            "limitations": [
                "Cluster bootstrap intervals concern synthetic map sampling",
                "Only 40 independently generated maps per profile",
                "Within-map coverage is estimated from 80 Poisson trials",
                "No inference about real HCT116 biological coverage",
            ],
        }, indent=2) + "\n"
    )

    print("\nSaved:", OUTPUT)
    print(
        "Bootstrap units are maps, not individual observations."
    )


if __name__ == "__main__":
    main()
