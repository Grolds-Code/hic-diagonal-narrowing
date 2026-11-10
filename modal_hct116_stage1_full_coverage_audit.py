"""Stage 1: spatial coverage audit of the complete fine scan."""

import hashlib
import io
from pathlib import Path

import modal

app = modal.App("hic-stage1-full-coverage-audit")

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
WINDOW = 64
MAX_DISTANCE = 64

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
def audit_windows(windows):
    import numpy as np

    path = Path(
        f"/research/matrices/{SOURCE}.npz"
    )

    raw = path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()

    if digest != EXPECTED_SHA:
        raise RuntimeError("Source checksum mismatch")

    with np.load(io.BytesIO(raw), allow_pickle=False) as z:
        x = z["bin_x_bp"].astype(np.int64)
        y = z["bin_y_bp"].astype(np.int64)
        counts = z["counts"].astype(float)

    d = (y - x) // RES
    selected = (d >= 1) & (d <= MAX_DISTANCE)

    if (
        len(x) != 45_831
        or int(selected.sum()) != 43_918
        or np.any(~np.isfinite(counts))
        or np.any(counts < 0)
        or np.any((y - x) % RES)
        or np.any(x > y)
    ):
        raise RuntimeError("Invalid source records")

    x = x[selected]
    y = y[selected]
    counts = counts[selected]

    bins = (END - START) // RES

    partners = np.zeros(bins, dtype=np.int64)
    signals = np.zeros(bins, dtype=float)

    for endpoint in (x, y):
        inside = (
            (endpoint >= START)
            & (endpoint < END)
        )

        idx = (
            (endpoint[inside] - START) // RES
        ).astype(int)

        np.add.at(partners, idx, 1)
        np.add.at(signals, idx, counts[inside])

    results = []

    for lo, hi in windows:
        if (
            lo < START + MAX_DISTANCE * RES
            or hi > END - MAX_DISTANCE * RES
            or hi - lo != WINDOW * RES
            or lo % RES
            or hi % RES
        ):
            raise RuntimeError(
                f"Invalid interior window: {lo}-{hi}"
            )

        start_bin = (lo - START) // RES
        end_bin = start_bin + WINDOW

        p = partners[start_bin:end_bin]
        s = signals[start_bin:end_bin]

        if (
            len(p) != WINDOW
            or np.any(p > 2 * MAX_DISTANCE)
        ):
            raise RuntimeError("Invalid endpoint support")

        zero = p == 0

        current = 0
        longest = 0

        for is_zero in zero:
            current = current + 1 if is_zero else 0
            longest = max(longest, current)

        results.append({
            "start_bp": lo,
            "end_bp": hi,
            "zero_contact_bins": int(zero.sum()),
            "longest_zero_run_bins": longest,
            "min_observed_partners": int(p.min()),
            "p10_observed_partners": float(
                np.percentile(p, 10)
            ),
            "median_observed_partners": float(
                np.median(p)
            ),
            "median_endpoint_contact_signal": float(
                np.median(s)
            ),
            "zero_contact_bin_starts_bp": ";".join(
                str(lo + int(i) * RES)
                for i in np.flatnonzero(zero)
            ),
        })

    return {
        "source_sha256": digest,
        "results": results,
    }


@app.local_entrypoint()
def main():
    import numpy as np
    import pandas as pd

    root = Path("research_artifacts/stage1")

    scored = pd.read_csv(
        root / "hct116_stage1_500bp_step_scored.csv"
    )

    previous = pd.read_csv(
        root / "hct116_stage1_endpoint_support_coordinates.csv"
    ).fillna("")

    regions = pd.read_csv(
        root / "hct116_stage1_candidate_deliverable_v1.csv"
    )

    if (
        len(scored) != 1937
        or int(scored["candidate_call"].sum()) != 93
        or len(previous) != 24
        or len(regions) != 5
    ):
        raise RuntimeError("Unexpected Stage 1 inputs")

    expected_starts = np.arange(
        START, END - WINDOW * RES + 1, RES
    )

    if (
        not np.array_equal(
            scored["start_bp"].to_numpy(),
            expected_starts,
        )
        or not np.array_equal(
            scored["end_bp"].to_numpy(),
            expected_starts + WINDOW * RES,
        )
    ):
        raise RuntimeError("Fine-scan geometry changed")

    assessable = (
        (scored["start_bp"] >= START + 32000)
        & (scored["end_bp"] <= END - 32000)
    )

    interior = scored.loc[assessable]

    if len(interior) != 1809:
        raise RuntimeError("Unexpected interior count")

    windows = [
        (int(r.start_bp), int(r.end_bp))
        for r in interior.itertuples(index=False)
    ]

    response = audit_windows.remote(windows)

    if response["source_sha256"] != EXPECTED_SHA:
        raise RuntimeError("Provenance mismatch")

    measured = pd.DataFrame(response["results"])

    if len(measured) != len(interior):
        raise RuntimeError("Incomplete coverage audit")

    result = scored.merge(
        measured,
        on=["start_bp", "end_bp"],
        how="left",
        validate="one_to_one",
        sort=False,
    )

    result["coverage_geometry"] = np.where(
        assessable.to_numpy(),
        "interior_complete",
        "boundary_limited",
    )

    if result.loc[
        assessable, "p10_observed_partners"
    ].isna().any():
        raise RuntimeError("Missing interior measurements")

    if result.loc[
        ~assessable, "p10_observed_partners"
    ].notna().any():
        raise RuntimeError("Boundary assessment not isolated")

    # Independent reproduction of earlier candidate and
    # control endpoint-support measurements.
    for old in previous.itertuples(index=False):
        matched = measured.loc[
            (measured["start_bp"] == int(old.start_bp))
            & (measured["end_bp"] == int(old.end_bp))
        ]

        if len(matched) != 1:
            raise RuntimeError(
                f"Previous window unavailable: {old.label}"
            )

        new = matched.iloc[0]

        for column in (
            "zero_contact_bins",
            "longest_zero_run_bins",
            "min_observed_partners",
            "p10_observed_partners",
            "median_observed_partners",
            "median_endpoint_contact_signal",
        ):
            if not np.isclose(
                float(new[column]),
                float(getattr(old, column)),
                rtol=0,
                atol=1e-9,
            ):
                raise RuntimeError(
                    f"{old.label}: {column} mismatch"
                )

    output = (
        root / "hct116_stage1_full_scan_coverage_v1.csv"
    )

    if output.exists():
        raise RuntimeError(
            f"Refusing to overwrite: {output}"
        )

    result.to_csv(output, index=False)

    print("\nSTAGE 1 — FULL-SCAN SPATIAL COVERAGE")
    print("=" * 80)
    print("Total scanned windows:", len(result))
    print("Complete interior assessments:", len(interior))
    print("Boundary-limited windows:", int((~assessable).sum()))
    print("Previous 24 support measurements reproduced: PASS")

    called = result.loc[result["candidate_call"]]

    if not called["coverage_geometry"].eq(
        "interior_complete"
    ).all():
        raise RuntimeError(
            "Candidate call near incomplete boundary"
        )

    print("\nCALLED WINDOWS — DESCRIPTIVE QC")
    print("-" * 80)

    for region in regions.itertuples(index=False):
        subset = called.loc[
            (
                called["start_bp"]
                >= region.footprint_start_bp
            )
            & (
                called["end_bp"]
                <= region.footprint_end_bp
            )
        ]

        print(
            f"{region.candidate_id}: "
            f"calls={len(subset)}, "
            f"windows_with_zero_bins="
            f"{int((subset['zero_contact_bins'] > 0).sum())}, "
            f"windows_with_zero_run_2plus="
            f"{int((subset['longest_zero_run_bins'] >= 2).sum())}, "
            f"minimum_p10_partners="
            f"{subset['p10_observed_partners'].min():.1f}"
        )

    print("\nGLOBAL INTERIOR COVERAGE DISTRIBUTION")
    print("-" * 80)

    for col in (
        "zero_contact_bins",
        "longest_zero_run_bins",
        "p10_observed_partners",
    ):
        values = result.loc[assessable, col]
        print(
            f"{col}: "
            f"min={values.min():.2f}, "
            f"median={values.median():.2f}, "
            f"max={values.max():.2f}"
        )

    print("\nLIMITATIONS")
    print("- No QC filtering or candidate reclassification.")
    print("- Overlapping windows are not independent samples.")
    print("- Contact support is measured within 32 kb.")
    print("- Zero contacts do not establish a deletion.")
    print("- Formal genomic QC remains unassessed.")
    print("Saved:", output)
    print("\nStep 375: PASS")
