"""Audit W80 estimands across independent fixed contact-intensity maps.

The analysis distinguishes:
1. A reference W80 formed from the generating mean profile.
2. The actual W80 of each fixed heterogeneous intensity map.

Independent Poisson observations are repeatedly generated
from each fixed map.

The independent unit for between-map variation is the map,
not each repeated Poisson observation.
"""

import json
from pathlib import Path

import numpy as np

from hic_diagonal_narrowing.poisson_uncertainty import (
    cutoff_exposure_w80_interval_from_bands,
)

SEED = 20261028
MAPS = 40
REPLICATES_PER_MAP = 80
K = 64
DMAX = 64

PROFILES = {
    "Region_A_like": (8.4, 2478),
    "Region_B_like": (7.5, 4695),
}


def width(profile):
    profile = np.asarray(profile, dtype=float)

    if np.any(profile < 0) or profile.sum() <= 0:
        raise ValueError("Invalid intensity profile")

    return int(
        np.searchsorted(
            np.cumsum(profile),
            0.8 * profile.sum(),
        ) + 1
    )


def run_profile(name, tau, expected_total, seed):
    rng = np.random.default_rng(seed)

    distances = np.arange(1, DMAX + 1)
    shape = np.exp(-distances / tau)

    # Expected count per genomic pair.
    mu = (
        shape * expected_total
        / (K * shape.sum())
    )

    population_w80 = width(K * mu)
    exposure = np.full(DMAX, K, dtype=np.int64)

    map_results = []

    for map_id in range(MAPS):
        # K independent pair-specific random intensities
        # within every separation band.
        #
        # Gamma shape 0.5, mean mu[d].
        # Once drawn, intensities are FIXED for all 80
        # observations from this particular map.
        pair_rates = rng.gamma(
            shape=0.5,
            scale=mu[None, :] / 0.5,
            size=(K, DMAX),
        )

        band_rates = pair_rates.sum(axis=0)
        map_w80 = width(band_rates)

        latent_hits = 0
        population_hits = 0
        intervals = []

        for _ in range(REPLICATES_PER_MAP):
            # Independent Poisson counts given the fixed map.
            counts = rng.poisson(
                band_rates
            ).astype(np.int64)

            result = (
                cutoff_exposure_w80_interval_from_bands(
                    counts,
                    exposure,
                    confidence_level=0.95,
                )
            )

            lo = int(result["cutoff_lower"])
            hi = int(result["cutoff_upper"])

            if not (
                1 <= lo
                <= int(result["observed_w80"])
                <= hi <= DMAX
            ):
                raise RuntimeError("Invalid confidence interval")

            latent_hits += int(
                lo <= map_w80 <= hi
            )

            population_hits += int(
                lo <= population_w80 <= hi
            )

            intervals.append(hi - lo)

        map_results.append({
            "map_id": map_id,
            "population_w80": population_w80,
            "fixed_map_w80": map_w80,
            "latent_target_shift_bins": (
                map_w80 - population_w80
            ),
            "conditional_latent_inclusion": (
                latent_hits / REPLICATES_PER_MAP
            ),
            "reference_population_inclusion": (
                population_hits / REPLICATES_PER_MAP
            ),
            "median_interval_span": float(
                np.median(intervals)
            ),
        })

    shifted = [
        row["latent_target_shift_bins"]
        for row in map_results
    ]

    latent_rates = [
        row["conditional_latent_inclusion"]
        for row in map_results
    ]

    population_rates = [
        row["reference_population_inclusion"]
        for row in map_results
    ]

    return {
        "region": name,
        "population_w80": population_w80,
        "independent_maps": MAPS,
        "poisson_observations_per_map": REPLICATES_PER_MAP,
        "total_observations": MAPS * REPLICATES_PER_MAP,
        "maps_with_shifted_latent_w80": int(
            np.count_nonzero(shifted)
        ),
        "fraction_maps_with_shifted_latent_w80": float(
            np.mean(np.asarray(shifted) != 0)
        ),
        "map_w80_range": [
            min(row["fixed_map_w80"] for row in map_results),
            max(row["fixed_map_w80"] for row in map_results),
        ],
        "mean_within_map_latent_inclusion": float(
            np.mean(latent_rates)
        ),
        "mean_reference_population_inclusion": float(
            np.mean(population_rates)
        ),
        "minimum_within_map_latent_inclusion": float(
            np.min(latent_rates)
        ),
        "median_within_map_latent_inclusion": float(
            np.median(latent_rates)
        ),
        "median_reference_population_inclusion": float(
            np.median(population_rates)
        ),
        "maps": map_results,
    }


def main():
    results = {}

    print("\nFIXED HETEROGENEOUS CONTACT-MAP AUDIT")
    print("=" * 65)

    print("Independent intensity maps:", MAPS)
    print("Poisson datasets per map:", REPLICATES_PER_MAP)
    print("Nominal confidence level: 95%")

    for index, (name, (tau, total)) in enumerate(
        PROFILES.items()
    ):
        result = run_profile(
            name, tau, total, SEED + index
        )

        results[name] = result

        print(f"\n{name}")
        print(
            "Reference population W80:",
            result["population_w80"],
        )
        print(
            "Fixed-map W80 range:",
            result["map_w80_range"],
        )
        print(
            "Fraction of maps with shifted target:",
            f"{result['fraction_maps_with_shifted_latent_w80']:.1%}",
        )
        print(
            "Mean within-map latent inclusion:",
            f"{result['mean_within_map_latent_inclusion']:.1%}",
        )
        print(
            "Mean reference population inclusion:",
            f"{result['mean_reference_population_inclusion']:.1%}",
        )
        print(
            "Minimum within-map latent inclusion:",
            f"{result['minimum_within_map_latent_inclusion']:.1%}",
        )

    output = Path(
        "research_artifacts/stage1/"
        "w80_fixed_heterogeneity_map_audit_v1.json"
    )

    if output.exists():
        raise RuntimeError(
            "Output already exists; refusing to overwrite"
        )

    output.parent.mkdir(parents=True, exist_ok=True)

    output.write_text(
        json.dumps({
            "experiment": "fixed_heterogeneity_map_audit",
            "seed": SEED,
            "maps_per_profile": MAPS,
            "replicates_per_map": REPLICATES_PER_MAP,
            "gamma_pair_shape": 0.5,
            "nominal_confidence": 0.95,
            "results": results,
            "limitations": (
                "Synthetic fixed heterogeneous intensities. "
                "Map-level variation is not biological replication. "
                "The reference population profile is a generating "
                "distribution expectation rather than the realized "
                "intensity profile of an individual map."
            ),
        }, indent=2) + "\n"
    )

    print("\nSaved:", output)
    print(
        "Interpretation: conditional Poisson sampling versus "
        "variation across independent synthetic intensity maps."
    )


if __name__ == "__main__":
    main()
