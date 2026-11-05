"""Real 500-bp HCT116 Hi-C contact maps with exploratory calls."""

import csv
import hashlib
import io
import json
from pathlib import Path

import modal

app = modal.App("hic-stage1-real-contact-map-calls")

volume = modal.Volume.from_name(
    "hic-diagonal-narrowing-data"
)

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install(
        "numpy==2.5.4",
        "matplotlib==3.10.7",
    )
    .add_local_file(
        "research_artifacts/stage1/hct116_stage1_merged_regions.csv",
        remote_path="/root/merged_regions.csv",
    )
)

SOURCE = "hct116_chr10_20000000_21000000_500bp_midpoint_32kb"

START = 20_000_000
END = 21_000_000
RES = 500
DMAX = 64


@app.function(
    image=image,
    volumes={"/research": volume},
    cpu=2,
    memory=4096,
    timeout=600,
)
def make_figures():
    import numpy as np

    import matplotlib
    matplotlib.use("Agg")

    import matplotlib.pyplot as plt

    from matplotlib.patches import Rectangle
    from matplotlib.lines import Line2D

    source_path = Path(
        f"/research/matrices/{SOURCE}.npz"
    )
    metadata_path = Path(
        f"/research/metadata/{SOURCE}.json"
    )

    raw = source_path.read_bytes()
    meta = json.loads(metadata_path.read_text())

    source_hash = hashlib.sha256(raw).hexdigest()

    if source_hash != meta["sha256"]:
        raise RuntimeError("Original contact source checksum mismatch")

    with np.load(io.BytesIO(raw), allow_pickle=False) as data:
        x = data["bin_x_bp"].astype(np.int64)
        y = data["bin_y_bp"].astype(np.int64)
        counts = data["counts"].astype(float)

    if not (len(x) == len(y) == len(counts) == 45831):
        raise RuntimeError("Unexpected source record count")

    if (
        np.any(x > y)
        or np.any(counts < 0)
        or np.any(~np.isfinite(counts))
    ):
        raise RuntimeError("Invalid contact records")

    d = (y - x) // RES

    offdiag = (d >= 1) & (d <= DMAX)

    if int(offdiag.sum()) != 43918:
        raise RuntimeError("Unexpected off-diagonal record count")

    regions = list(csv.DictReader(
        Path("/root/merged_regions.csv").open()
    ))

    if len(regions) != 2:
        raise RuntimeError("Expected two merged regions")

    expected = [
        (20_296_000, 20_336_000),
        (20_720_000, 20_792_000),
    ]

    for row, (left, right) in zip(regions, expected):
        if (
            int(row["start_bp"]) != left
            or int(row["end_bp"]) != right
            or row["qc_status"] != "unassessed"
        ):
            raise RuntimeError("Candidate region mismatch")

    cmap = plt.get_cmap("magma").copy()
    cmap.set_bad("#e9edf2")

    vmax = float(np.log1p(
        np.percentile(counts[offdiag], 99)
    ))

    def contact_matrix(lo, hi):
        if (
            lo < START or hi > END
            or (hi - lo) % RES
        ):
            raise RuntimeError("Invalid requested map extent")

        n = (hi - lo) // RES
        matrix = np.full((n, n), np.nan)

        # The diagonal is excluded from the W80 analysis.
        # Distances greater than 32 kb were not extracted.
        for distance in range(1, min(DMAX, n - 1) + 1):
            rows = np.arange(n - distance)
            cols = rows + distance

            mid = lo + (
                rows + distance / 2
            ) * RES

            valid = (mid >= START) & (mid < END)

            matrix[rows[valid], cols[valid]] = 0.0
            matrix[cols[valid], rows[valid]] = 0.0

        selected = (
            offdiag
            & (x >= lo)
            & (y < hi)
        )

        ii = ((x[selected] - lo) // RES).astype(int)
        jj = ((y[selected] - lo) // RES).astype(int)
        cc = counts[selected]

        if np.any(~np.isfinite(matrix[ii, jj])):
            raise RuntimeError("Contact outside eligible map support")

        np.add.at(matrix, (ii, jj), cc)
        np.add.at(matrix, (jj, ii), cc)

        if not np.allclose(
            matrix,
            matrix.T,
            equal_nan=True,
        ):
            raise RuntimeError("Contact map is not symmetric")

        return np.ma.masked_invalid(np.log1p(matrix))

    def add_boxes(ax, region, labels=False):
        left = int(region["start_bp"]) / 1e6
        right = int(region["end_bp"]) / 1e6

        peak_left = int(region["peak_start_bp"]) / 1e6
        peak_right = int(region["peak_end_bp"]) / 1e6

        ax.add_patch(Rectangle(
            (left, left),
            right - left,
            right - left,
            fill=False,
            edgecolor="#00bcd4",
            linewidth=2.3,
            zorder=5,
        ))

        ax.add_patch(Rectangle(
            (peak_left, peak_left),
            peak_right - peak_left,
            peak_right - peak_left,
            fill=False,
            edgecolor="#a5ff69",
            linewidth=2,
            linestyle="--",
            zorder=6,
        ))

    def show_map(ax, lo, hi, label):
        data = contact_matrix(lo, hi)

        im = ax.imshow(
            data,
            origin="lower",
            extent=(
                lo / 1e6, hi / 1e6,
                lo / 1e6, hi / 1e6,
            ),
            cmap=cmap,
            vmin=0,
            vmax=vmax,
            interpolation="nearest",
            aspect="equal",
        )

        ax.set_title(label, fontsize=11)
        ax.set_xlabel("chr10 position (Mb)")
        ax.set_ylabel("chr10 position (Mb)")
        ax.tick_params(labelsize=8)

        return im

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "pdf.fonttype": 42,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

    outputs = {}

    # Figure 1: full 1-Mb contact-map overview.
    fig, ax = plt.subplots(figsize=(8.4, 7.7))

    im = show_map(
        ax,
        START,
        END,
        "HCT116 chr10 | 500-bp Hi-C contact map",
    )

    for region in regions:
        add_boxes(ax, region)

    fig.colorbar(
        im,
        ax=ax,
        fraction=0.046,
        pad=0.035,
        label="log(1 + raw contact count)",
    )

    fig.text(
        0.5, 0.015,
        "500-bp bins | 32-kb detector windows | 8-kb scan step | "
        "off-diagonal contacts up to 32 kb",
        ha="center",
        fontsize=8,
    )

    fig.tight_layout(rect=(0, 0.035, 1, 1))

    for extension in ("png", "pdf"):
        stream = io.BytesIO()
        fig.savefig(
            stream,
            format=extension,
            dpi=230,
            bbox_inches="tight",
        )
        outputs[
            f"hct116_chr10_500bp_contact_map_overview.{extension}"
        ] = stream.getvalue()

    plt.close(fig)

    # Figure 2: zoomed real contact maps.
    zooms = [
        (20_216_000, 20_408_000),
        (20_656_000, 20_848_000),
    ]

    fig, axes = plt.subplots(
        1, 2,
        figsize=(12.6, 6.2),
    )

    for idx, (ax, region, bounds) in enumerate(
        zip(axes, regions, zooms)
    ):
        im = show_map(
            ax,
            bounds[0],
            bounds[1],
            (
                f"Region {'AB'[idx]} | "
                f"peak W80 = {float(region['peak_w80_bins']) * 0.5:g} kb | "
                f"score = {float(region['peak_score']):.3f}"
            ),
        )

        add_boxes(ax, region)

    fig.colorbar(
        im,
        ax=axes,
        fraction=0.025,
        pad=0.03,
        label="log(1 + raw contact count)",
    )

    legend = [
        Line2D(
            [0], [0],
            color="#00bcd4",
            lw=2.5,
            label="Merged detection footprint",
        ),
        Line2D(
            [0], [0],
            color="#71ae25",
            lw=2.5,
            linestyle="--",
            label="Peak detection window",
        ),
    ]

    fig.legend(
        handles=legend,
        loc="lower center",
        ncol=2,
        bbox_to_anchor=(0.48, 0.04),
        frameon=False,
        fontsize=10,
    )

    fig.suptitle(
        "Observed 500-bp HCT116 Hi-C maps with exploratory "
        "diagonal-narrowing calls",
        fontsize=13,
        y=0.97,
    )

    fig.text(
        0.48, 0.025,
        "Boxes locate detector windows, not exact chromatin "
        "boundaries. QC unassessed.",
        ha="center",
        fontsize=8,
    )

    fig.subplots_adjust(
        left=0.07,
        right=0.87,
        bottom=0.20,
        top=0.87,
        wspace=0.30,
    )

    for extension in ("png", "pdf"):
        stream = io.BytesIO()
        fig.savefig(
            stream,
            format=extension,
            dpi=250,
            bbox_inches="tight",
        )
        outputs[
            f"hct116_chr10_500bp_candidate_contact_maps.{extension}"
        ] = stream.getvalue()

    plt.close(fig)

    return {
        "source_sha256": source_hash,
        "offdiagonal_records": int(offdiag.sum()),
        "merged_regions": len(regions),
        "files": outputs,
    }


@app.local_entrypoint()
def main():
    result = make_figures.remote()

    output = Path("figures/stage1")
    output.mkdir(parents=True, exist_ok=True)

    for filename, content in result["files"].items():
        path = output / filename

        if path.exists():
            raise RuntimeError(
                f"Refusing to overwrite: {path}"
            )

        path.write_bytes(content)

        print(f"Saved: {path} ({len(content):,} bytes)")

    print("\nREAL HCT116 CONTACT MAP VERIFICATION")
    print("=" * 60)
    print("Source SHA-256:", result["source_sha256"])
    print("Off-diagonal records:", result["offdiagonal_records"])
    print("Merged regions:", result["merged_regions"])
    print("Bin resolution: 500 bp")
    print("Detection window: 32 kb")
    print("Scan step: 8 kb")
    print("Biological interpretation: EXPLORATORY")
    print("Step 349: PASS")
