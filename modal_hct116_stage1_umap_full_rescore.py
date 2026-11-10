"""Full fine-grid mappability sensitivity after coverage triage."""

import hashlib
import io
import subprocess
from pathlib import Path

import modal

app = modal.App("hic-stage1-umap-full-rescore")

volume = modal.Volume.from_name(
    "hic-diagonal-narrowing-data"
)

image = (
    modal.Image.debian_slim(python_version="3.12")
    .apt_install("ca-certificates", "curl")
    .pip_install("numpy==2.5.4")
    .run_commands(
        "curl -fL --retry 3 --max-time 120 "
        "-o /usr/local/bin/bigWigToBedGraph "
        "https://hgdownload.soe.ucsc.edu/admin/exe/"
        "linux.x86_64/bigWigToBedGraph",
        "chmod +x /usr/local/bin/bigWigToBedGraph",
    )
)

START = 20_000_000
END = 21_000_000
RES = 500
N = 2000
W = 64
DMAX = 64

SOURCE = (
    "hct116_chr10_20000000_21000000_"
    "500bp_midpoint_32kb"
)

SHA = (
    "fb01bc848b4fecb342a450b8d87165393"
    "a89de019fc27c82e72d8519c07f5de5"
)

ROOT = (
    "https://hgdownload.soe.ucsc.edu/gbdb/"
    "hg38/hoffmanMappability/"
)

TRACKS = {
    "k50": ROOT + "k50.Umap.MultiTrackMappability.bw",
    "k100": ROOT + "k100.Umap.MultiTrackMappability.bw",
}


@app.function(
    image=image,
    volumes={"/research": volume},
    cpu=2,
    memory=4096,
    timeout=900,
)
def compute():
    import numpy as np

    source = Path(
        f"/research/matrices/{SOURCE}.npz"
    )
    raw = source.read_bytes()

    if hashlib.sha256(raw).hexdigest() != SHA:
        raise RuntimeError("Contact source checksum mismatch")

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

    x = x[selected]
    y = y[selected]
    counts = counts[selected]
    distance = distance[selected]

    # Interior scan windows only. Their contact endpoints
    # remain inside the extracted one-megabase interval.
    starts = np.arange(64, N - W - 64 + 1)
    ends = starts + W

    if len(starts) != 1809:
        raise RuntimeError("Interior geometry mismatch")

    in_bounds = (x >= START) & (y < END)
    x = x[in_bounds]
    y = y[in_bounds]
    counts = counts[in_bounds]
    distance = distance[in_bounds]

    xi = (x - START) // RES
    yi = (y - START) // RES
    mid = (xi + yi) // 2

    def evaluate(valid_bins):
        signal = np.zeros((N, DMAX), dtype=np.int64)
        eligible = np.zeros((N, DMAX), dtype=np.int64)

        valid_observed = (
            valid_bins[xi] & valid_bins[yi]
        )

        np.add.at(
            signal,
            (
                mid[valid_observed],
                distance[valid_observed] - 1,
            ),
            counts[valid_observed],
        )

        positions = np.arange(N)

        for d in range(1, DMAX + 1):
            left = positions - d // 2
            right = left + d

            inside = (
                (left >= 0)
                & (right < N)
            )

            eligible[inside, d - 1] = (
                valid_bins[left[inside]]
                & valid_bins[right[inside]]
            )

        signal_prefix = np.vstack([
            np.zeros((1, DMAX), dtype=np.int64),
            np.cumsum(signal, axis=0),
        ])

        eligible_prefix = np.vstack([
            np.zeros((1, DMAX), dtype=np.int64),
            np.cumsum(eligible, axis=0),
        ])

        band_signal = (
            signal_prefix[ends]
            - signal_prefix[starts]
        )

        band_eligible = (
            eligible_prefix[ends]
            - eligible_prefix[starts]
        )

        totals = band_signal.sum(axis=1)
        valid_pairs = band_eligible.sum(axis=1)
        minimum_pairs = band_eligible.min(axis=1)

        good = (
            (minimum_pairs > 0)
            & (totals > 0)
        )

        width = np.full(len(starts), -1, dtype=int)

        means = (
            band_signal[good]
            / band_eligible[good]
        )

        cumulative = np.cumsum(means, axis=1)
        targets = 0.8 * cumulative[:, -1]

        width[good] = (
            np.argmax(
                cumulative >= targets[:, None],
                axis=1,
            ) + 1
        )

        invalid_prefix = np.r_[
            0, np.cumsum(~valid_bins)
        ]

        masked_bins = (
            invalid_prefix[ends]
            - invalid_prefix[starts]
        )

        return {
            "width": width.tolist(),
            "totals": totals.tolist(),
            "valid_pairs": valid_pairs.tolist(),
            "minimum_pairs": minimum_pairs.tolist(),
            "masked_bins": masked_bins.tolist(),
        }

    scenarios = {
        "unmasked": evaluate(
            np.ones(N, dtype=bool)
        )
    }

    for name, url in TRACKS.items():
        process = subprocess.run(
            [
                "/usr/local/bin/bigWigToBedGraph",
                "-chrom=chr10",
                f"-start={START}",
                f"-end={END}",
                url,
                "stdout",
            ],
            capture_output=True,
            text=True,
            timeout=240,
            check=False,
        )

        if process.returncode != 0:
            raise RuntimeError(
                f"{name}: {process.stderr[-500:]}"
            )

        base_values = np.zeros(
            END - START, dtype=np.float32
        )

        interval_count = 0
        previous_end = START

        for line in process.stdout.splitlines():
            fields = line.split()

            if len(fields) != 4 or fields[0] != "chr10":
                raise RuntimeError(
                    f"{name}: invalid BedGraph record"
                )

            left = int(fields[1])
            right = int(fields[2])
            value = float(fields[3])

            if (
                left < previous_end
                or right <= left
                or not np.isfinite(value)
                or not 0 <= value <= 1
            ):
                raise RuntimeError(
                    f"{name}: invalid Umap interval"
                )

            previous_end = right

            a = max(left, START) - START
            b = min(right, END) - START

            if b > a:
                base_values[a:b] = value
                interval_count += 1

        if interval_count == 0:
            raise RuntimeError(
                f"{name}: no valid track intervals"
            )

        means = base_values.reshape(N, RES).mean(
            axis=1
        )

        for threshold in (0.50, 0.80, 0.95):
            key = f"{name}_{int(threshold * 100):02d}"

            result = evaluate(means >= threshold)
            result["genomic_masked_bins"] = int(
                (means < threshold).sum()
            )

            scenarios[key] = result

    return {
        "starts_bp": (
            START + starts * RES
        ).tolist(),
        "scenarios": scenarios,
        "source_sha256": SHA,
    }


@app.local_entrypoint()
def main():
    import numpy as np
    import pandas as pd

    from hic_diagonal_narrowing.detect import (
        detect_narrowing
    )

    root = Path("research_artifacts/stage1")

    profile = pd.read_csv(
        root / "hct116_stage1_500bp_step_profile.csv"
    )

    triage = pd.read_csv(
        root / "hct116_stage1_coverage_triage_v1.csv"
    )

    regions = pd.read_csv(
        root / "hct116_stage1_candidate_deliverable_v1.csv"
    )

    eligibility_audit = Path(
        "research_artifacts/stage2/"
        "hct116_stage2_umap_eligibility_audit_v1.csv"
    )

    destination = (
        root / "hct116_stage1_umap_full_rescore_v1.csv"
    )

    if destination.exists():
        raise RuntimeError(
            f"Refusing to overwrite: {destination}"
        )

    if (
        len(profile) != 1937
        or len(triage) != 1937
        or len(regions) != 5
    ):
        raise RuntimeError("Unexpected Stage 1 inputs")

    if not profile[
        ["start_bp", "end_bp"]
    ].equals(
        triage[["start_bp", "end_bp"]]
    ):
        raise RuntimeError("Stage 1 rows do not align")

    response = compute.remote()

    if response["source_sha256"] != SHA:
        raise RuntimeError("Source provenance mismatch")

    interior = (
        (profile["start_bp"] >= START + 32000)
        & (profile["end_bp"] <= END - 32000)
    ).to_numpy()

    expected_starts = profile.loc[
        interior, "start_bp"
    ].astype(int).tolist()

    if expected_starts != response["starts_bp"]:
        raise RuntimeError(
            "Remote/local window geometry mismatch"
        )

    baseline = response["scenarios"]["unmasked"]

    if not np.array_equal(
        np.asarray(baseline["width"]),
        profile.loc[interior, "w80_bins"].to_numpy(
            dtype=int
        ),
    ):
        raise RuntimeError("Unmasked W80 mismatch")

    if not np.array_equal(
        np.asarray(baseline["totals"]),
        profile.loc[interior, "total_contacts"].to_numpy(
            dtype=int
        ),
    ):
        raise RuntimeError("Unmasked count mismatch")

    if not all(
        value == 4096
        for value in baseline["valid_pairs"]
    ):
        raise RuntimeError("Unmasked eligibility mismatch")

    print("\nSTAGE 1 — FULL UMAP RESCORING")
    print("=" * 90)
    print("Original HCT116 data reproduced: PASS")
    print("Interior windows:", len(expected_starts))

    audit = (
        pd.read_csv(eligibility_audit)
        if eligibility_audit.exists()
        else None
    )

    outputs = []

    for name, result in response["scenarios"].items():
        work = profile.copy()

        width = np.asarray(
            result["width"], dtype=float
        )
        width[width < 0] = np.nan

        totals = np.asarray(
            result["totals"], dtype=float
        )

        valid_pairs = np.asarray(
            result["valid_pairs"], dtype=int
        )

        minimum_pairs = np.asarray(
            result["minimum_pairs"], dtype=int
        )

        work.loc[interior, "w80_bins"] = width
        work.loc[interior, "w80_bp"] = width * RES
        work.loc[interior, "total_contacts"] = totals

        valid_masked = (
            np.isfinite(width)
            & (minimum_pairs > 0)
            & (totals > 0)
        )

        allowed = (
            triage["coverage_flag"]
            .eq("no_coverage_warning")
            .to_numpy()
        )

        allowed = allowed.copy()
        allowed[interior] &= valid_masked

        work["qc_status"] = np.where(
            allowed, "unassessed", "coverage_warning"
        )

        scored, footprints = detect_narrowing(
            work,
            score_threshold=0.20,
            min_reference_windows=3,
        )

        current = scored.loc[interior].copy()

        current["scenario"] = name
        current["valid_pairs"] = valid_pairs
        current["minimum_pairs"] = minimum_pairs
        current["masked_bins"] = result["masked_bins"]
        current["removed_pair_fraction"] = (
            1 - valid_pairs / 4096
        )
        current["formal_qc_status"] = "unassessed"

        outputs.append(current)

        if name == "unmasked":
            if (
                int(scored["candidate_call"].sum()) != 25
                or len(footprints) != 2
            ):
                raise RuntimeError(
                    "Stage 383 baseline not reproduced"
                )

        if audit is not None and name != "unmasked":
            track, threshold_text = name.split("_")
            threshold = int(threshold_text) / 100

            for label in ("R1", "R2"):
                match = audit.loc[
                    (audit["candidate_id"] == label)
                    & (audit["role"] == "candidate")
                    & (audit["track"] == track)
                    & np.isclose(
                        audit["threshold"], threshold
                    )
                ]

                if len(match) != 1:
                    raise RuntimeError(
                        "Stage 2 eligibility row missing"
                    )

                item = match.iloc[0]
                candidate_start = int(
                    item["window_center_bp"] - 16000
                )

                verified = current.loc[
                    current["start_bp"] == candidate_start
                ]

                if len(verified) != 1:
                    raise RuntimeError(
                        "Candidate window missing"
                    )

                if int(
                    verified.iloc[0]["valid_pairs"]
                ) != int(item["valid_pairs"]):
                    raise RuntimeError(
                        f"{name} {label}: "
                        "Umap eligibility disagreement"
                    )

        called = int(scored["candidate_call"].sum())

        print(f"\n{name}")
        print("-" * 65)
        print("Eligible windows:", int(allowed.sum()))
        print("Candidate calls:", called)
        print("Merged footprints:", len(footprints))

        for region in regions.itertuples(index=False):
            member = (
                (scored["start_bp"] >=
                 region.footprint_start_bp)
                & (scored["end_bp"] <=
                   region.footprint_end_bp)
            )

            retained = scored.loc[
                member & scored["candidate_call"]
            ]

            print(
                f"  {region.candidate_id}: "
                f"calls_in_original_footprint="
                f"{len(retained)}"
            )

    full = pd.concat(
        outputs, ignore_index=True
    )

    if len(full) != 1809 * 7:
        raise RuntimeError("Incomplete scenario results")

    full.to_csv(destination, index=False)

    print("\nVERIFICATION")
    print("Unmasked source measurements: PASS")
    print("Coverage-screened baseline: PASS")
    print(
        "Existing Umap candidate eligibility cross-check:",
        "PASS" if audit is not None else "NOT AVAILABLE"
    )
    print("Full masked reference recalculation: COMPLETE")
    print("Read-length assumptions: k50 and k100")
    print("Umap thresholds: SENSITIVITY SETTINGS")
    print("Formal QC: UNASSESSED")
    print("Saved:", destination)
    print("\nStep 386: PASS")
