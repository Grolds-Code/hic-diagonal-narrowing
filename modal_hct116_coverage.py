"""Read-only coverage diagnostics for the HCT116 Hi-C pilot matrix."""

import hashlib
import json
from pathlib import Path

import modal

app = modal.App("hic-diagonal-narrowing-coverage")

volume = modal.Volume.from_name(
    "hic-diagonal-narrowing-data"
)

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install("numpy==2.5.4")
)

BASENAME = "hct116_chr10_20000000_20100000_500bp"


@app.function(
    image=image,
    volumes={"/research": volume},
    cpu=1,
    memory=1024,
    timeout=120,
)
def analyze_coverage():
    import numpy as np

    matrix_path = Path(
        f"/research/matrices/{BASENAME}.npz"
    )
    metadata_path = Path(
        f"/research/metadata/{BASENAME}.json"
    )

    metadata = json.loads(metadata_path.read_text())

    digest = hashlib.sha256(
        matrix_path.read_bytes()
    ).hexdigest()

    if digest != metadata["sha256"]:
        raise RuntimeError("Matrix checksum mismatch")

    with np.load(matrix_path, allow_pickle=False) as data:
        matrix = data["contacts"]

    if matrix.shape != (200, 200):
        raise RuntimeError("Unexpected matrix shape")

    if not np.array_equal(matrix, matrix.T):
        raise RuntimeError("Matrix is not symmetric")

    if not np.all(np.isfinite(matrix)):
        raise RuntimeError("Nonfinite contact counts")

    if np.any(matrix < 0):
        raise RuntimeError("Negative contact counts")

    unique_total = float(np.triu(matrix).sum())

    if not np.isclose(
        unique_total,
        metadata["unique_contact_signal"]
    ):
        raise RuntimeError("Contact total mismatch")

    row_totals = matrix.sum(axis=1)

    # Each band includes one copy of every contact pair.
    bands = [
        (0, 1),
        (1, 5),
        (5, 10),
        (10, 20),
        (20, 50),
        (50, 100),
        (100, 200),
    ]

    band_results = []

    for low, high in bands:
        possible = 0
        observed = 0
        counts = 0.0

        for offset in range(low, high):
            diagonal = np.diagonal(
                matrix,
                offset=offset
            )

            possible += diagonal.size
            observed += int(
                np.count_nonzero(diagonal)
            )
            counts += float(diagonal.sum())

        band_results.append({
            "separation_bins": (
                str(low) if high == low + 1
                else f"{low}-{high - 1}"
            ),
            "possible_pairs": possible,
            "observed_pairs": observed,
            "occupancy": observed / possible,
            "contact_signal": counts,
        })

    if not np.isclose(
        sum(b["contact_signal"] for b in band_results),
        unique_total
    ):
        raise RuntimeError("Band totals do not reconcile")

    return {
        "checksum_verified": True,
        "shape": list(matrix.shape),
        "resolution_bp": metadata["resolution_bp"],
        "unique_contact_signal": unique_total,
        "zero_signal_bins": int(
            np.count_nonzero(row_totals == 0)
        ),
        "row_signal_percentiles": [
            float(x) for x in np.percentile(
                row_totals,
                [0, 10, 25, 50, 75, 90, 100]
            )
        ],
        "bands": band_results,
    }


@app.local_entrypoint()
def main():
    result = analyze_coverage.remote()

    print("\nHCT116 CONTACT COVERAGE DIAGNOSTICS")
    print("=" * 65)

    print("Checksum:", "PASS" if result[
        "checksum_verified"
    ] else "FAIL")
    print("Shape:", result["shape"])
    print("Resolution:", result["resolution_bp"], "bp")
    print(
        "Unique contact signal:",
        result["unique_contact_signal"]
    )
    print(
        "Zero-signal bins:",
        result["zero_signal_bins"]
    )

    print("\nPER-BIN SIGNAL PERCENTILES")
    print("Percentiles: 0, 10, 25, 50, 75, 90, 100")
    print("Values:", result["row_signal_percentiles"])

    print("\nDISTANCE-DEPENDENT OCCUPANCY")
    print(
        f"{'Bins':<13}"
        f"{'Possible':>10}"
        f"{'Observed':>11}"
        f"{'Occupancy':>13}"
        f"{'Signal':>12}"
    )

    for band in result["bands"]:
        print(
            f"{band['separation_bins']:<13}"
            f"{band['possible_pairs']:>10}"
            f"{band['observed_pairs']:>11}"
            f"{band['occupancy']:>13.4f}"
            f"{band['contact_signal']:>12.1f}"
        )

    print("\nRead-only coverage diagnostics: PASS")
