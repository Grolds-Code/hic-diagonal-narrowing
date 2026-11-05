"""Test sensitivity of HCT116 candidate W80 to zero-contact bins."""

import hashlib
import io
from pathlib import Path

import modal

app = modal.App("hic-stage1-zero-bin-sensitivity")

volume = modal.Volume.from_name(
    "hic-diagonal-narrowing-data"
)

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install("numpy==2.5.4")
)

START = 20_000_000
RES = 500
WINDOW = 64
DMAX = 64

SOURCE = (
    "hct116_chr10_20000000_21000000_"
    "500bp_midpoint_32kb"
)

EXPECTED_SHA = (
    "fb01bc848b4fecb342a450b8d87165393"
    "a89de019fc27c82e72d8519c07f5de5"
)


@app.function(
    image=image,
    volumes={"/research": volume},
    cpu=2,
    memory=4096,
    timeout=600,
)
def evaluate(windows):
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

    separation = (y - x) // RES
    selected = (separation >= 1) & (separation <= DMAX)

    if len(x) != 45831 or int(selected.sum()) != 43918:
        raise RuntimeError("Source record count mismatch")

    records = {
        (int(a), int(b)): int(c)
        for a, b, c in zip(
            x[selected], y[selected], counts[selected]
        )
    }

    if len(records) != 43918:
        raise RuntimeError("Duplicate contact coordinates")

    def calculate_w80(bands):
        total = bands.sum()

        if total <= 0:
            raise RuntimeError("No contact signal")

        return int(
            np.searchsorted(
                np.cumsum(bands),
                0.8 * total,
            ) + 1
        )

    results = []

    for window in windows:
        lo = int(window["start_bp"])
        hi = int(window["end_bp"])
        zero = set(window["zero_bin_starts"])

        if hi - lo != WINDOW * RES:
            raise RuntimeError("Unexpected window width")

        if any(
            z < lo or z >= hi or z % RES
            for z in zero
        ):
            raise RuntimeError("Invalid zero-bin coordinate")

        start_bin = (lo - START) // RES

        original_bands = []
        masked_bands = []
        remaining_per_band = []

        observed_pairs = 0
        total_contacts = 0
        masked_pairs = 0

        for d in range(1, DMAX + 1):
            # Exactly 64 pairs have midpoints inside a
            # 64-bin window at each separation.
            first = start_bin - d // 2
            left_bins = np.arange(
                first, first + WINDOW
            )

            left_bp = START + RES * left_bins
            right_bp = left_bp + RES * d

            values = np.array([
                records.get((int(a), int(b)), 0)
                for a, b in zip(left_bp, right_bp)
            ], dtype=np.int64)

            observed_pairs += int(
                np.count_nonzero(values)
            )
            total_contacts += int(values.sum())

            # Mask only the reported zero-contact bins.
            valid = np.array([
                int(a) not in zero and int(b) not in zero
                for a, b in zip(left_bp, right_bp)
            ])

            n_valid = int(valid.sum())

            if n_valid == 0:
                raise RuntimeError(
                    "A distance band lost all eligible pairs"
                )

            if np.any(values[~valid] != 0):
                raise RuntimeError(
                    "Reported zero-contact bin has signal"
                )

            masked_pairs += WINDOW - n_valid
            remaining_per_band.append(n_valid)

            original_bands.append(
                float(values.sum()) / WINDOW
            )
            masked_bands.append(
                float(values[valid].sum()) / n_valid
            )

        if observed_pairs != window["observed_pairs"]:
            raise RuntimeError(
                f"{window['label']}: observed-pair mismatch"
            )

        if total_contacts != window["total_contacts"]:
            raise RuntimeError(
                f"{window['label']}: contact-total mismatch"
            )

        original_w80 = calculate_w80(
            np.array(original_bands)
        )

        if original_w80 != window["w80_bins"]:
            raise RuntimeError(
                f"{window['label']}: original W80 mismatch"
            )

        masked_w80 = calculate_w80(
            np.array(masked_bands)
        )

        results.append({
            "region": window["label"],
            "start_bp": lo,
            "end_bp": hi,
            "zero_contact_bins": len(zero),
            "original_w80_bp": original_w80 * RES,
            "masked_w80_bp": masked_w80 * RES,
            "w80_change_bp": (
                masked_w80 - original_w80
            ) * RES,
            "masked_eligible_pairs": masked_pairs,
            "remaining_eligible_pairs": (
                WINDOW * DMAX - masked_pairs
            ),
            "min_pairs_per_distance": min(
                remaining_per_band
            ),
            "original_total_contacts": total_contacts,
            "original_observed_pairs": observed_pairs,
            "formal_qc_status": "unassessed",
        })

    return results


@app.local_entrypoint()
def main():
    import numpy as np
    import pandas as pd

    root = Path("research_artifacts/stage1")

    calls = pd.read_csv(
        root / "hct116_stage1_500bp_step_calls.csv"
    )
    regions = pd.read_csv(
        root / "hct116_stage1_500bp_step_regions.csv"
    )
    support = pd.read_csv(
        root / "hct116_stage1_endpoint_support_coordinates.csv"
    ).fillna("")

    if (len(calls), len(regions)) != (93, 5):
        raise RuntimeError("Unexpected candidate inputs")

    windows = []

    for i, region in enumerate(regions.itertuples(), 1):
        label = f"R{i}"

        group = calls.loc[
            (calls["start_bp"] >= region.start_bp)
            & (calls["end_bp"] <= region.end_bp)
        ]

        best = group["narrowing_score"].max()

        tied = group.loc[
            np.isclose(
                group["narrowing_score"].to_numpy(),
                best,
                rtol=0,
                atol=1e-12,
            )
        ].sort_values("start_bp")

        representative = tied.iloc[len(tied) // 2]

        row = support.loc[
            (support["label"] == label)
            & (support["kind"] == "candidate")
        ]

        if len(row) != 1:
            raise RuntimeError(
                f"Missing support record: {label}"
            )

        row = row.iloc[0]

        if (
            int(row["start_bp"])
            != int(representative["start_bp"])
        ):
            raise RuntimeError(
                f"Representative window mismatch: {label}"
            )

        zero_text = str(row["zero_contact_bin_starts_bp"])

        zero_bins = (
            [int(v) for v in zero_text.split(";")]
            if zero_text else []
        )

        if len(zero_bins) != int(row["zero_contact_bins"]):
            raise RuntimeError(
                f"Zero-bin count mismatch: {label}"
            )

        windows.append({
            "label": label,
            "start_bp": int(representative["start_bp"]),
            "end_bp": int(representative["end_bp"]),
            "w80_bins": int(representative["w80_bins"]),
            "total_contacts": int(
                representative["total_contacts"]
            ),
            "observed_pairs": int(
                representative["observed_pairs"]
            ),
            "zero_bin_starts": zero_bins,
        })

    result = pd.DataFrame(evaluate.remote(windows))

    if len(result) != 5:
        raise RuntimeError("Incomplete sensitivity analysis")

    destination = (
        root / "hct116_stage1_zero_bin_sensitivity.csv"
    )

    if destination.exists():
        raise RuntimeError(
            f"Refusing to overwrite: {destination}"
        )

    result.to_csv(destination, index=False)

    print("\nSTAGE 1 — ZERO-BIN SENSITIVITY")
    print("=" * 90)

    print(
        result[
            [
                "region",
                "zero_contact_bins",
                "original_w80_bp",
                "masked_w80_bp",
                "w80_change_bp",
                "masked_eligible_pairs",
                "min_pairs_per_distance",
            ]
        ].to_string(index=False)
    )

    print("\nVERIFICATION")
    print("Original contact totals and W80: PASS")
    print("Distance-specific eligible pairs: RECALCULATED")
    print("Missing-bin classification: NOT ESTABLISHED")
    print("Formal QC status: UNASSESSED")
    print("Saved:", destination)
    print("Step 361: PASS")
