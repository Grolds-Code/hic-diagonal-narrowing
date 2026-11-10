"""Audit low-but-nonzero endpoint coverage in Stage 1."""

import hashlib
import io
from pathlib import Path

import modal

app = modal.App("hic-stage1-low-support-audit")

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

SHA = (
    "fb01bc848b4fecb342a450b8d87165393"
    "a89de019fc27c82e72d8519c07f5de5"
)

START = 20_000_000
END = 21_000_000
RES = 500


@app.function(
    image=image,
    volumes={"/research": volume},
    cpu=2,
    memory=4096,
    timeout=600,
)
def evaluate(windows):
    import numpy as np

    path = Path(
        f"/research/matrices/{SOURCE}.npz"
    )

    raw = path.read_bytes()

    if hashlib.sha256(raw).hexdigest() != SHA:
        raise RuntimeError("HCT116 checksum mismatch")

    with np.load(io.BytesIO(raw), allow_pickle=False) as z:
        x = z["bin_x_bp"].astype(np.int64)
        y = z["bin_y_bp"].astype(np.int64)
        counts = z["counts"]

    d = (y - x) // RES
    selected = (d >= 1) & (d <= 64)

    if (
        len(x) != 45831
        or int(selected.sum()) != 43918
        or np.any(x > y)
        or np.any((y - x) % RES)
        or np.any(~np.isfinite(counts))
        or np.any(counts < 0)
    ):
        raise RuntimeError("Invalid source records")

    x = x[selected]
    y = y[selected]

    partners = np.zeros(
        (END - START) // RES, dtype=np.int64
    )

    for endpoints in (x, y):
        inside = (
            (endpoints >= START)
            & (endpoints < END)
        )

        indices = (
            (endpoints[inside] - START) // RES
        ).astype(int)

        np.add.at(partners, indices, 1)

    def longest_run(mask):
        current = 0
        best = 0

        for flag in mask:
            current = current + 1 if flag else 0
            best = max(best, current)

        return best

    records = []

    for start, end in windows:
        start = int(start)
        end = int(end)

        if (
            start < START + 32000
            or end > END - 32000
            or end - start != 32000
            or start % RES
        ):
            raise RuntimeError("Invalid interior window")

        first = (start - START) // RES

        # The eligible midpoint contact pairs can have
        # endpoints 16 kb outside either window boundary.
        central = partners[first:first + 64]
        endpoint = partners[first - 32:first + 96]

        if len(central) != 64 or len(endpoint) != 128:
            raise RuntimeError("Invalid endpoint geometry")

        median = float(np.median(endpoint))

        if median <= 0:
            raise RuntimeError("No local partner support")

        zero = endpoint == 0
        low10 = endpoint < 0.10 * median
        low25 = endpoint < 0.25 * median

        records.append({
            "start_bp": start,
            "end_bp": end,
            "central_zero_bins_check": int(
                (central == 0).sum()
            ),
            "central_p10_check": float(
                np.percentile(central, 10)
            ),
            "endpoint_median_partners": median,
            "endpoint_p10_partners": float(
                np.percentile(endpoint, 10)
            ),
            "endpoint_p10_to_median_ratio": float(
                np.percentile(endpoint, 10) / median
            ),
            "endpoint_zero_bins": int(zero.sum()),
            "endpoint_longest_zero_run": longest_run(
                zero
            ),
            "endpoint_low10_bins": int(low10.sum()),
            "endpoint_low10_longest_run": longest_run(
                low10
            ),
            "endpoint_low25_bins": int(low25.sum()),
            "endpoint_low25_longest_run": longest_run(
                low25
            ),
        })

    return records


@app.local_entrypoint()
def main():
    import numpy as np
    import pandas as pd

    root = Path("research_artifacts/stage1")

    audited = pd.read_csv(
        root / "hct116_stage1_full_scan_coverage_v1.csv"
    )

    regions = pd.read_csv(
        root / "hct116_stage1_candidate_deliverable_v1.csv"
    )

    interior = audited.loc[
        audited["coverage_geometry"].eq(
            "interior_complete"
        )
    ].copy()

    if (
        len(audited) != 1937
        or len(interior) != 1809
        or len(regions) != 5
    ):
        raise RuntimeError("Unexpected Stage 1 inputs")

    windows = list(zip(
        interior["start_bp"].astype(int),
        interior["end_bp"].astype(int),
    ))

    measured = pd.DataFrame(
        evaluate.remote(windows)
    )

    result = interior.merge(
        measured,
        on=["start_bp", "end_bp"],
        how="inner",
        validate="one_to_one",
    )

    if len(result) != 1809:
        raise RuntimeError("Incomplete audit")

    if not np.array_equal(
        result["central_zero_bins_check"].to_numpy(),
        result["zero_contact_bins"].to_numpy()
    ):
        raise RuntimeError(
            "Previous zero-bin measurements disagree"
        )

    if not np.allclose(
        result["central_p10_check"],
        result["p10_observed_partners"],
        rtol=0,
        atol=1e-10,
    ):
        raise RuntimeError(
            "Previous partner-support values disagree"
        )

    output = (
        root /
        "hct116_stage1_low_support_full_scan_v1.csv"
    )

    if output.exists():
        raise RuntimeError(
            f"Refusing to overwrite: {output}"
        )

    result.to_csv(output, index=False)

    print("\nSTAGE 1 — LOCAL ENDPOINT COVERAGE AUDIT")
    print("=" * 82)
    print("Audited interior windows:", len(result))
    print("Original coverage measurements reproduced: PASS")
    print("Endpoint span: 64 kb")
    print("Low-support settings: 10% and 25% of local median")

    for region in regions.itertuples(index=False):
        group = result.loc[
            result["candidate_call"].astype(bool)
            & (
                result["start_bp"]
                >= region.footprint_start_bp
            )
            & (
                result["end_bp"]
                <= region.footprint_end_bp
            )
        ]

        if len(group) != int(region.supporting_windows):
            raise RuntimeError(
                f"{region.candidate_id}: call mismatch"
            )

        print(f"\n{region.candidate_id}")
        print("-" * 65)
        print("Original calls:", len(group))
        print(
            "Endpoint zero runs >= 2 bins:",
            int(
                (
                    group["endpoint_longest_zero_run"] >= 2
                ).sum()
            )
        )
        print(
            "Endpoint low10 runs >= 2 bins:",
            int(
                (
                    group["endpoint_low10_longest_run"] >= 2
                ).sum()
            )
        )
        print(
            "Endpoint low25 runs >= 2 bins:",
            int(
                (
                    group["endpoint_low25_longest_run"] >= 2
                ).sum()
            )
        )
        print(
            "Minimum endpoint p10/median ratio:",
            round(
                group[
                    "endpoint_p10_to_median_ratio"
                ].min(),
                4,
            )
        )

    noncalls = result.loc[
        ~result["candidate_call"].astype(bool)
    ]

    print("\nNON-CALLED INTERIOR WINDOWS")
    print("-" * 65)
    print("Windows:", len(noncalls))

    for column in (
        "endpoint_longest_zero_run",
        "endpoint_low10_longest_run",
        "endpoint_low25_longest_run",
    ):
        print(
            f"{column}: "
            f"{int((noncalls[column] >= 2).sum())} "
            "windows with runs >= 2"
        )

    print("\nLIMITATIONS")
    print("- Low-support cutoffs are exploratory.")
    print("- Overlapping windows are not independent.")
    print("- Local coverage may reflect biological copy number.")
    print("- Missing-bin identity is not established.")
    print("- No detector calls have been changed.")
    print("- Formal QC remains unassessed.")
    print("Saved:", output)
    print("\nStep 380: PASS")
