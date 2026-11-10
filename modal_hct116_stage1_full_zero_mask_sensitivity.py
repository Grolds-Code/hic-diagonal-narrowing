"""Test zero-endpoint masking across all 93 Stage 1 calls."""

import csv
import hashlib
import io
from pathlib import Path

import modal

app = modal.App("hic-stage1-full-zero-mask-sensitivity")

volume = modal.Volume.from_name(
    "hic-diagonal-narrowing-data"
)

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install("numpy==2.5.4")
)

SOURCE = (
    "hct116_chr10_20000000_21000000_"
    "500bp_midpoint_32kb"
)

EXPECTED_SHA = (
    "fb01bc848b4fecb342a450b8d87165393"
    "a89de019fc27c82e72d8519c07f5de5"
)

RES = 500
WINDOW = 64
DMAX = 64


@app.function(
    image=image,
    volumes={"/research": volume},
    cpu=2,
    memory=4096,
    timeout=600,
)
def evaluate(windows, zero_bins):
    import numpy as np

    source = Path(
        f"/research/matrices/{SOURCE}.npz"
    )

    raw = source.read_bytes()

    if hashlib.sha256(raw).hexdigest() != EXPECTED_SHA:
        raise RuntimeError("Source checksum mismatch")

    with np.load(io.BytesIO(raw), allow_pickle=False) as z:
        x = z["bin_x_bp"].astype(np.int64)
        y = z["bin_y_bp"].astype(np.int64)
        counts = z["counts"].astype(np.int64)

    distance = (y - x) // RES
    selected = (distance >= 1) & (distance <= DMAX)

    if (
        len(x) != 45831
        or int(selected.sum()) != 43918
        or np.any(x > y)
        or np.any((y - x) % RES)
        or np.any(counts < 0)
    ):
        raise RuntimeError("Unexpected source geometry")

    lookup = {
        (int(a), int(b)): int(c)
        for a, b, c in zip(
            x[selected], y[selected], counts[selected]
        )
    }

    if len(lookup) != 43918:
        raise RuntimeError("Duplicate contact records")

    zero = np.array(sorted(zero_bins), dtype=np.int64)
    result = []

    for window in windows:
        start = int(window["start_bp"])
        end = int(window["end_bp"])

        if end - start != WINDOW * RES:
            raise RuntimeError("Unexpected window width")

        original_bands = []
        masked_means = []

        removed_pairs = 0
        minimum_valid = WINDOW
        observed_pairs = 0

        for d in range(1, DMAX + 1):
            left = (
                start
                - (d // 2) * RES
                + np.arange(WINDOW, dtype=np.int64) * RES
            )
            right = left + d * RES

            values = np.fromiter(
                (
                    lookup.get((int(a), int(b)), 0)
                    for a, b in zip(left, right)
                ),
                dtype=np.int64,
                count=WINDOW,
            )

            valid = (
                ~np.isin(left, zero)
                & ~np.isin(right, zero)
            )

            if np.any(values[~valid] != 0):
                raise RuntimeError(
                    "Zero-support bin has observed contact"
                )

            n_valid = int(valid.sum())

            if n_valid == 0:
                raise RuntimeError(
                    "Distance lost all eligible pairs"
                )

            removed_pairs += WINDOW - n_valid
            minimum_valid = min(minimum_valid, n_valid)
            observed_pairs += int(np.count_nonzero(values))

            original_bands.append(int(values.sum()))
            masked_means.append(
                float(values[valid].sum()) / n_valid
            )

        original_bands = np.asarray(
            original_bands, dtype=np.int64
        )
        masked_means = np.asarray(
            masked_means, dtype=float
        )

        total = int(original_bands.sum())

        original_w80 = (
            int(np.searchsorted(
                np.cumsum(original_bands),
                0.8 * total,
            )) + 1
        ) * RES

        masked_total = masked_means.sum()

        if masked_total <= 0:
            raise RuntimeError("No masked contact signal")

        masked_w80 = (
            int(np.searchsorted(
                np.cumsum(masked_means),
                0.8 * masked_total,
            )) + 1
        ) * RES

        if (
            original_w80 != int(window["w80_bp"])
            or total != int(window["total_contacts"])
            or observed_pairs != int(window["observed_pairs"])
        ):
            raise RuntimeError(
                f"Baseline mismatch at {start}"
            )

        result.append({
            "region": window["region"],
            "start_bp": start,
            "end_bp": end,
            "original_w80_bp": original_w80,
            "zero_masked_w80_bp": masked_w80,
            "w80_change_bp": masked_w80 - original_w80,
            "pairs_touching_zero_bins": removed_pairs,
            "affected_pair_fraction": (
                removed_pairs / (WINDOW * DMAX)
            ),
            "minimum_valid_pairs_per_distance": minimum_valid,
            "formal_qc_status": "unassessed",
        })

    return result


@app.local_entrypoint()
def main():
    import pandas as pd

    root = Path("research_artifacts/stage1")

    scored = pd.read_csv(
        root / "hct116_stage1_500bp_step_scored.csv"
    )

    regions = pd.read_csv(
        root / "hct116_stage1_candidate_deliverable_v1.csv"
    )

    support = pd.read_csv(
        root / "hct116_stage1_full_scan_coverage_v1.csv"
    )

    if (
        len(scored) != 1937
        or len(regions) != 5
        or len(support) != 1937
    ):
        raise RuntimeError("Unexpected Stage 1 inputs")

    calls = scored.loc[
        scored["candidate_call"].astype(str)
        .str.lower().eq("true")
    ].copy()

    if len(calls) != 93:
        raise RuntimeError("Expected 93 original calls")

    zero_positions = set()

    for entry in support[
        "zero_contact_bin_starts_bp"
    ].dropna():
        for token in str(entry).split(";"):
            if token.strip():
                zero_positions.add(int(token))

    if len(zero_positions) != 8:
        raise RuntimeError(
            "Unexpected number of zero-support bins"
        )

    windows = []

    for call in calls.itertuples(index=False):
        matching = regions.loc[
            (regions["footprint_start_bp"] <= call.start_bp)
            & (regions["footprint_end_bp"] >= call.end_bp)
        ]

        if len(matching) != 1:
            raise RuntimeError("Ambiguous candidate membership")

        windows.append({
            "region": matching.iloc[0]["candidate_id"],
            "start_bp": int(call.start_bp),
            "end_bp": int(call.end_bp),
            "w80_bp": int(call.w80_bp),
            "total_contacts": int(call.total_contacts),
            "observed_pairs": int(call.observed_pairs),
        })

    result = pd.DataFrame(
        evaluate.remote(windows, sorted(zero_positions))
    )

    if len(result) != 93:
        raise RuntimeError("Incomplete sensitivity analysis")

    output = (
        root /
        "hct116_stage1_called_zero_endpoint_mask_sensitivity_v1.csv"
    )

    if output.exists():
        raise RuntimeError(
            f"Refusing to overwrite: {output}"
        )

    result.to_csv(output, index=False)

    print("\nSTAGE 1 — FULL CALLED-WINDOW MASK SENSITIVITY")
    print("=" * 82)
    print("Original candidate windows:", len(result))
    print("Identified zero-support bins:", len(zero_positions))
    print("Original W80 and contact totals reproduced: PASS")

    for region, group in result.groupby(
        "region", sort=True
    ):
        changed = group.loc[
            group["w80_change_bp"] != 0
        ]

        print(f"\n{region}")
        print("-" * 60)
        print("Original calls:", len(group))
        print(
            "Windows exposed to zero bins:",
            int(
                (group["pairs_touching_zero_bins"] > 0).sum()
            )
        )
        print(
            "Windows with changed W80:",
            len(changed)
        )
        print(
            "W80 change range:",
            int(group["w80_change_bp"].min()),
            "to",
            int(group["w80_change_bp"].max()),
            "bp"
        )
        print(
            "Affected pair fraction:",
            f"{group['affected_pair_fraction'].min():.2%}",
            "to",
            f"{group['affected_pair_fraction'].max():.2%}"
        )

    print("\nLIMITATIONS")
    print("- This is a sensitivity test, not a QC filter.")
    print("- Only observed zero-support bins were excluded.")
    print("- Low-but-nonzero coverage was not masked.")
    print("- Mappability-based filtering is a separate test.")
    print("- Candidate scores have not been recalculated.")
    print("- No biological deletion has been established.")
    print("- Formal QC remains unassessed.")
    print("Saved:", output)
    print("\nStep 378: PASS")
