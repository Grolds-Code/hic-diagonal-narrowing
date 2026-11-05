"""Descriptive per-bin contact-support audit for Stage 1."""

import hashlib
import io
from pathlib import Path

import modal

app = modal.App("hic-stage1-endpoint-support")

volume = modal.Volume.from_name(
    "hic-diagonal-narrowing-data"
)

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install("numpy==2.5.4")
)

START = 20_000_000
END = 21_000_000
RES = 500
DMAX = 64
WINDOW = 32_000

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
def audit(windows):
    import numpy as np

    source = Path(
        f"/research/matrices/{SOURCE}.npz"
    )

    raw = source.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()

    if digest != EXPECTED_SHA:
        raise RuntimeError("Unexpected source checksum")

    with np.load(io.BytesIO(raw), allow_pickle=False) as data:
        x = data["bin_x_bp"].astype(np.int64)
        y = data["bin_y_bp"].astype(np.int64)
        counts = data["counts"].astype(float)

    distances = (y - x) // RES
    selected = (distances >= 1) & (distances <= DMAX)

    if (
        len(x) != 45_831
        or int(selected.sum()) != 43_918
        or np.any(~np.isfinite(counts))
        or np.any(counts < 0)
        or np.any(x > y)
        or np.any((y - x) % RES != 0)
    ):
        raise RuntimeError("Invalid source contact records")

    x = x[selected]
    y = y[selected]
    counts = counts[selected]

    n = (END - START) // RES

    observed_partners = np.zeros(n, dtype=np.int64)
    endpoint_signal = np.zeros(n, dtype=float)

    # Count observed off-diagonal contact partners and
    # raw signal for each genomic bin.
    for endpoint in (x, y):
        inside = (
            (endpoint >= START)
            & (endpoint < END)
        )

        indices = (
            (endpoint[inside] - START) // RES
        ).astype(int)

        np.add.at(observed_partners, indices, 1)
        np.add.at(endpoint_signal, indices, counts[inside])

    # Every window used below is at least DMAX bins
    # from the target-region boundary. Thus each bin
    # has all 128 potential off-diagonal partners
    # within the extracted 32-kb separation range.
    if np.any(observed_partners > 2 * DMAX):
        raise RuntimeError(
            "Unexpected duplicate or excessive partner counts"
        )

    results = []

    for label, kind, lo, hi in windows:
        if (
            lo < START + DMAX * RES
            or hi > END - DMAX * RES
            or hi - lo != WINDOW
            or lo % RES
            or hi % RES
        ):
            raise RuntimeError(
                f"Invalid window geometry: {label}"
            )

        first = (lo - START) // RES
        last = (hi - START) // RES

        partners = observed_partners[first:last]
        signals = endpoint_signal[first:last]

        if len(partners) != 64:
            raise RuntimeError("Incorrect window bin count")

        zero = partners == 0

        longest_run = 0
        current_run = 0

        for missing in zero:
            if missing:
                current_run += 1
                longest_run = max(
                    longest_run, current_run
                )
            else:
                current_run = 0

        results.append({
            "label": label,
            "kind": kind,
            "start_bp": lo,
            "end_bp": hi,
            "zero_contact_bins": int(zero.sum()),
            "zero_contact_bin_starts_bp": ";".join(
                str(lo + int(j) * RES)
                for j in np.flatnonzero(zero)
            ),
            "longest_zero_run_bins": longest_run,
            "min_observed_partners": int(partners.min()),
            "p10_observed_partners": float(
                np.percentile(partners, 10)
            ),
            "median_observed_partners": float(
                np.median(partners)
            ),
            "median_endpoint_contact_signal": float(
                np.median(signals)
            ),
            "missing_bin_status": "not_assessed",
        })

    return results


@app.local_entrypoint()
def main():
    import numpy as np
    import pandas as pd

    root = Path("research_artifacts/stage1")

    baseline = pd.read_csv(
        root / "w80_profile_v2_offdiag.csv"
    )
    calls = pd.read_csv(
        root / "hct116_stage1_500bp_step_calls.csv"
    )
    regions = pd.read_csv(
        root / "hct116_stage1_500bp_step_regions.csv"
    )

    if (len(baseline), len(calls), len(regions)) != (
        122, 93, 5
    ):
        raise RuntimeError("Unexpected source profiles")

    windows = []

    # Representative tied-peak window for each region.
    for i, region in enumerate(regions.itertuples(), 1):
        group = calls.loc[
            (calls["start_bp"] >= region.start_bp)
            & (calls["end_bp"] <= region.end_bp)
        ]

        best = group["narrowing_score"].max()

        tied = group.loc[
            np.isclose(
                group["narrowing_score"],
                best,
                atol=1e-12,
                rtol=0,
            )
        ].sort_values("start_bp")

        representative = tied.iloc[len(tied) // 2]

        windows.append((
            f"R{i}",
            "candidate",
            int(representative["start_bp"]),
            int(representative["end_bp"]),
        ))

    # Non-overlapping reference windows outside
    # all five candidate footprints.
    eligible = np.ones(len(baseline), dtype=bool)

    for region in regions.itertuples():
        eligible &= (
            (baseline["end_bp"].to_numpy() <= region.start_bp)
            | (baseline["start_bp"].to_numpy() >= region.end_bp)
        )

    pool = baseline.loc[eligible].sort_values("start_bp")
    last_end = -1
    control_count = 0

    for row in pool.itertuples():
        if (
            row.start_bp >= last_end
            and row.start_bp >= START + DMAX * RES
            and row.end_bp <= END - DMAX * RES
        ):
            control_count += 1
            windows.append((
                f"C{control_count}",
                "control",
                int(row.start_bp),
                int(row.end_bp),
            ))
            last_end = row.end_bp

    if len(windows) != 5 + control_count or control_count < 5:
        raise RuntimeError("Insufficient control windows")

    result = pd.DataFrame(audit.remote(windows))

    if len(result) != len(windows):
        raise RuntimeError("Incomplete support audit")

    candidates = result.loc[
        result["kind"] == "candidate"
    ]
    controls = result.loc[
        result["kind"] == "control"
    ]

    output = root / "hct116_stage1_endpoint_support_coordinates.csv"

    if output.exists():
        raise RuntimeError(
            f"Refusing to overwrite: {output}"
        )

    result.to_csv(output, index=False)

    print("\nSTAGE 1 — PER-BIN CONTACT SUPPORT")
    print("=" * 85)
    print("Original contact checksum: PASS")
    print("Candidate windows:", len(candidates))
    print("Non-overlapping interior controls:", len(controls))
    print("Bins per window: 64")
    print("Possible off-diagonal partners per bin: 128")

    print("\nCANDIDATE MEASUREMENTS")
    print(candidates[
        [
            "label",
            "start_bp",
            "zero_contact_bins",
            "longest_zero_run_bins",
            "min_observed_partners",
            "p10_observed_partners",
            "median_observed_partners",
            "median_endpoint_contact_signal",
        ]
    ].to_string(index=False))

    print("\nCONTROL DISTRIBUTION")
    for column in (
        "zero_contact_bins",
        "longest_zero_run_bins",
        "min_observed_partners",
        "median_observed_partners",
    ):
        print(
            f"{column}: "
            f"min={controls[column].min():.2f}, "
            f"median={controls[column].median():.2f}, "
            f"max={controls[column].max():.2f}"
        )

    print("\nLIMITATIONS")
    print("- Endpoint support is descriptive, not formal QC.")
    print("- Zero observed contacts do not prove a missing bin.")
    print("- Controls are not coverage-matched.")
    print("- Biological validity remains unassessed.")
    print("\nSaved:", output)
    print("Step 360: PASS")
