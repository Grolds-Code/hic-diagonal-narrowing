"""Plot five fine-grid calls on real 500-bp HCT116 Hi-C maps."""

import csv
import hashlib
import io
import json
from pathlib import Path

import modal

app = modal.App("hic-stage1-five-region-figures")

volume = modal.Volume.from_name(
    "hic-diagonal-narrowing-data"
)

ROOT = "research_artifacts/stage1"

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install("numpy==2.5.4", "matplotlib==3.10.7")
    .add_local_file(
        f"{ROOT}/hct116_stage1_500bp_step_profile.csv",
        remote_path="/root/profile.csv",
    )
    .add_local_file(
        f"{ROOT}/hct116_stage1_500bp_step_regions.csv",
        remote_path="/root/regions.csv",
    )
    .add_local_file(
        f"{ROOT}/hct116_stage1_fine_step_quality_summary.csv",
        remote_path="/root/quality.csv",
    )
)

START = 20_000_000
END = 21_000_000
RES = 500
DMAX = 64

SOURCE = (
    "hct116_chr10_20000000_21000000_"
    "500bp_midpoint_32kb"
)


@app.function(
    image=image,
    volumes={"/research": volume},
    cpu=2,
    memory=4096,
    timeout=600,
)
def render_figures():
    import numpy as np

    import matplotlib
    matplotlib.use("Agg")

    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    def read_csv(path):
        with open(path, newline="") as f:
            return list(csv.DictReader(f))

    profile = read_csv("/root/profile.csv")
    regions = read_csv("/root/regions.csv")
    quality = read_csv("/root/quality.csv")

    if (
        len(profile) != 1937
        or len(regions) != 5
        or len(quality) != 5
    ):
        raise RuntimeError("Unexpected fine-grid inputs")

    starts = np.array(
        [int(r["start_bp"]) for r in profile]
    )
    ends = np.array(
        [int(r["end_bp"]) for r in profile]
    )
    widths = np.array(
        [float(r["w80_bp"]) / 1000 for r in profile]
    )

    if (
        starts[0] != START
        or starts[-1] != END - 32_000
        or np.any(np.diff(starts) != RES)
        or np.any(ends - starts != 32_000)
        or not np.all(np.isfinite(widths))
    ):
        raise RuntimeError("Invalid fine-grid profile")

    if (
        sum(int(r["supporting_windows"]) for r in regions)
        != 93
    ):
        raise RuntimeError("Candidate count mismatch")

    for region, qc in zip(regions, quality):
        if (
            region["qc_status"] != "unassessed"
            or qc["qc_status"] != "unassessed"
            or int(region["start_bp"]) != int(qc["start_bp"])
            or int(region["end_bp"]) != int(qc["end_bp"])
        ):
            raise RuntimeError("Region/quality mismatch")

    source_path = Path(
        f"/research/matrices/{SOURCE}.npz"
    )
    meta_path = Path(
        f"/research/metadata/{SOURCE}.json"
    )

    raw = source_path.read_bytes()
    meta = json.loads(meta_path.read_text())
    sha = hashlib.sha256(raw).hexdigest()

    expected_sha = (
        "fb01bc848b4fecb342a450b8d87165393"
        "a89de019fc27c82e72d8519c07f5de5"
    )

    if sha != expected_sha or sha != meta["sha256"]:
        raise RuntimeError("Original Hi-C source mismatch")

    with np.load(io.BytesIO(raw), allow_pickle=False) as z:
        x = z["bin_x_bp"].astype(np.int64)
        y = z["bin_y_bp"].astype(np.int64)
        counts = z["counts"].astype(float)

    distances = (y - x) // RES
    offdiag = (distances >= 1) & (distances <= DMAX)

    if (
        len(x) != 45831
        or int(offdiag.sum()) != 43918
        or np.any(~np.isfinite(counts))
        or np.any(counts < 0)
    ):
        raise RuntimeError("Original contacts failed validation")

    vmax = float(
        np.log1p(np.percentile(counts[offdiag], 99))
    )

    cmap = plt.get_cmap("magma").copy()
    cmap.set_bad("#e5e9ef")

    def contact_map(lo, hi):
        if (
            lo < START or hi > END
            or lo % RES or hi % RES
            or lo >= hi
        ):
            raise ValueError("Invalid contact-map extent")

        n = (hi - lo) // RES
        matrix = np.full((n, n), np.nan, dtype=float)

        # Only off-diagonal separations <= 32 kb
        # are represented in the extracted source.
        for d in range(1, min(DMAX, n - 1) + 1):
            indices = np.arange(n - d)
            matrix[indices, indices + d] = 0
            matrix[indices + d, indices] = 0

        selected = (
            offdiag
            & (x >= lo)
            & (y < hi)
        )

        ii = ((x[selected] - lo) // RES).astype(int)
        jj = ((y[selected] - lo) // RES).astype(int)

        if np.any(~np.isfinite(matrix[ii, jj])):
            raise RuntimeError("Contact outside extracted support")

        np.add.at(matrix, (ii, jj), counts[selected])
        np.add.at(matrix, (jj, ii), counts[selected])

        if not np.allclose(
            matrix, matrix.T, equal_nan=True
        ):
            raise RuntimeError("Map symmetry failure")

        return np.ma.masked_invalid(np.log1p(matrix))

    def draw_map(ax, lo, hi):
        im = ax.imshow(
            contact_map(lo, hi),
            origin="lower",
            extent=[lo / 1e6, hi / 1e6,
                    lo / 1e6, hi / 1e6],
            cmap=cmap,
            vmin=0,
            vmax=vmax,
            interpolation="nearest",
            aspect="equal",
        )
        ax.set_xlabel("chr10 position (Mb)")
        ax.set_ylabel("chr10 position (Mb)")
        ax.tick_params(labelsize=8)
        return im

    def draw_footprint(ax, region):
        left = int(region["start_bp"]) / 1e6
        right = int(region["end_bp"]) / 1e6

        ax.add_patch(Rectangle(
            (left, left),
            right - left,
            right - left,
            fill=False,
            edgecolor="#00c7d9",
            linewidth=2.0,
            zorder=5,
        ))

    def save_figure(fig, basename, output):
        for extension in ("png", "pdf"):
            buffer = io.BytesIO()

            fig.savefig(
                buffer,
                format=extension,
                dpi=230,
                bbox_inches="tight",
                facecolor="white",
            )

            output[f"{basename}.{extension}"] = (
                buffer.getvalue()
            )

        plt.close(fig)

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "pdf.fonttype": 42,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

    outputs = {}

    # Figure 1: real 1-Mb Hi-C map and aligned W80 track.
    fig, (ax_map, ax_width) = plt.subplots(
        2, 1,
        figsize=(10, 12),
        gridspec_kw={"height_ratios": [4.5, 1]},
        layout="constrained",
    )

    im = draw_map(ax_map, START, END)

    for i, region in enumerate(regions, 1):
        draw_footprint(ax_map, region)

        center = (
            int(region["start_bp"])
            + int(region["end_bp"])
        ) / 2e6

        ax_map.text(
            center, center, f"R{i}",
            color="white",
            fontsize=8,
            fontweight="bold",
            ha="center",
            va="center",
            bbox={
                "facecolor": "#17202a",
                "edgecolor": "none",
                "alpha": 0.85,
                "pad": 2,
            },
            zorder=10,
        )

    fig.colorbar(
        im,
        ax=ax_map,
        shrink=0.80,
        pad=0.02,
        label="log(1 + raw contact count)",
    )

    centers = (starts + ends) / 2e6

    ax_width.plot(
        centers,
        widths,
        color="#234c78",
        linewidth=1.2,
        label="W80 from 500-bp-step profile",
    )

    for i, region in enumerate(regions, 1):
        left = int(region["start_bp"]) / 1e6
        right = int(region["end_bp"]) / 1e6

        ax_width.axvspan(
            left, right,
            color="#00aebd",
            alpha=0.18,
        )

        ax_width.text(
            (left + right) / 2,
            float(np.max(widths)) + 0.3,
            f"R{i}",
            ha="center",
            va="bottom",
            fontsize=8,
        )

    ax_width.set_xlim(START / 1e6, END / 1e6)
    ax_width.set_ylim(0, float(np.max(widths)) + 1.2)
    ax_width.set_xlabel("chr10 position (Mb)")
    ax_width.set_ylabel("W80 (kb)")
    ax_width.grid(alpha=0.2)
    ax_width.legend(loc="upper right", fontsize=8)

    fig.suptitle(
        "HCT116 | 500-bp Hi-C diagonal and "
        "fine-grid narrowing candidates",
        fontsize=14,
    )

    fig.text(
        0.5, 0.005,
        "32-kb measurement windows | 500-bp scan step | "
        "exploratory calls | QC unassessed",
        ha="center",
        fontsize=8,
    )

    save_figure(
        fig,
        "hct116_stage1_fine_scan_overview",
        outputs,
    )

    # Figure 2: five matched-scale candidate zooms.
    fig, axes = plt.subplots(
        2, 3,
        figsize=(15, 10),
        layout="constrained",
    )

    for i, (region, ax) in enumerate(
        zip(regions, axes.flat), 1
    ):
        left = int(region["start_bp"])
        right = int(region["end_bp"])

        center = (left + right) // 2
        lo = ((center - 60_000) // RES) * RES
        hi = lo + 120_000

        im = draw_map(ax, lo, hi)
        draw_footprint(ax, region)

        w80_kb = float(region["peak_w80_bins"]) * 0.5
        score = float(region["peak_score"])

        ax.set_title(
            f"R{i}: chr10:{left:,}-{right:,}\n"
            f"minimum W80 = {w80_kb:g} kb | "
            f"score = {score:.3f}",
            fontsize=10,
        )

    # Sixth panel explains what the figures establish.
    notes = axes.flat[5]
    notes.axis("off")

    notes.text(
        0.02, 0.95,
        "Reading these figures\n\n"
        "Cyan squares:\n"
        "Merged detector footprints\n\n"
        "Heatmap:\n"
        "Real, unnormalized HCT116 contacts\n"
        "within the 32-kb extracted band\n\n"
        "Gray regions:\n"
        "Unrepresented contact separations\n\n"
        "Input resolution: 500 bp\n"
        "Measurement window: 32 kb\n"
        "Scan step: 500 bp\n\n"
        "These are exploratory detections.\n"
        "They are not precise boundaries.\n"
        "Formal genomic QC is pending.",
        transform=notes.transAxes,
        va="top",
        fontsize=10,
        linespacing=1.4,
    )

    fig.colorbar(
        im,
        ax=list(axes.flat[:5]),
        location="bottom",
        shrink=0.65,
        pad=0.07,
        label="log(1 + raw contact count)",
    )

    fig.suptitle(
        "Five fine-grid HCT116 diagonal-narrowing candidates",
        fontsize=14,
    )

    save_figure(
        fig,
        "hct116_stage1_fine_scan_five_regions",
        outputs,
    )

    return {
        "sha256": sha,
        "windows": len(profile),
        "regions": len(regions),
        "files": outputs,
    }


@app.local_entrypoint()
def main():
    output = Path("figures/stage1")
    output.mkdir(parents=True, exist_ok=True)

    expected = [
        "hct116_stage1_fine_scan_overview.png",
        "hct116_stage1_fine_scan_overview.pdf",
        "hct116_stage1_fine_scan_five_regions.png",
        "hct116_stage1_fine_scan_five_regions.pdf",
    ]

    for name in expected:
        if (output / name).exists():
            raise RuntimeError(
                f"Refusing to overwrite: {output / name}"
            )

    result = render_figures.remote()

    if (
        result["windows"] != 1937
        or result["regions"] != 5
        or set(result["files"]) != set(expected)
    ):
        raise RuntimeError("Unexpected remote result")

    for name, content in result["files"].items():
        path = output / name
        path.write_bytes(content)
        print(f"Saved: {path} ({len(content):,} bytes)")

    print("\nSTAGE 1 — FIVE-REGION FIGURES")
    print("=" * 60)
    print("Original HCT116 SHA-256:", result["sha256"])
    print("Fine-grid windows:", result["windows"])
    print("Candidate regions:", result["regions"])
    print("Figure files:", len(result["files"]))
    print("Contact-map bin resolution: 500 bp")
    print("Measurement window: 32 kb")
    print("Scan step: 500 bp")
    print("Formal QC: UNASSESSED")
    print("Step 356: PASS")
