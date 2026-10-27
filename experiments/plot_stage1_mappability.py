"""Visualize HCT116 W80 sensitivity to Umap-based masking."""

import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

base = Path("research_artifacts/stage1")
source = base / "umap_mask_w80_v1.json"
meta = json.loads(
    (base / "umap_mask_w80_v1.metadata.json").read_text()
)
data = json.loads(source.read_text())

assert hashlib.sha256(
    source.read_bytes()
).hexdigest() == meta["output_sha256"]

assert data["assembly"] == "GRCh38"
assert data["resolution_bp"] == 500

regions = ["Region A", "Region B"]
colors = {
    "Region A": "#48668A",
    "Region B": "#B17D48",
}

baseline = {
    item["region"]: item["matched_narrowing_score"]
    for item in data["baseline"]
}

positions = np.array([0, 1, 2, 3.5, 4.5, 5.5])

labels = []
scores = {name: [] for name in regions}
pair_loss = {name: [] for name in regions}

for track_name in ["k50", "k100"]:
    scenarios = data["tracks"][track_name]["scenarios"]
    scenarios = sorted(
        scenarios,
        key=lambda item: item["validity_threshold"]
    )

    assert len(scenarios) == 3

    for scenario in scenarios:
        threshold = scenario["validity_threshold"]
        labels.append(f"{track_name}\n{threshold:.2f}")

        comparisons = {
            item["region"]: item
            for item in scenario["comparisons"]
        }

        assert set(comparisons) == set(regions)

        for name in regions:
            result = comparisons[name]

            scores[name].append(
                result["matched_narrowing_score"]
            )
            pair_loss[name].append(
                100 * result["peak_removed_pair_fraction"]
            )

assert len(labels) == 6

fig, (ax1, ax2) = plt.subplots(
    1, 2,
    figsize=(11.4, 5.9),
)

for name in regions:
    color = colors[name]

    ax1.axhline(
        baseline[name],
        color=color,
        linestyle=":",
        linewidth=1.4,
        alpha=0.8,
    )

    # Keep distinct Umap tracks visually disconnected.
    for sl in [slice(0, 3), slice(3, 6)]:
        ax1.plot(
            positions[sl],
            np.asarray(scores[name])[sl],
            marker="o",
            linewidth=2,
            markersize=5,
            color=color,
            label=name if sl.start == 0 else None,
        )

        ax2.plot(
            positions[sl],
            np.asarray(pair_loss[name])[sl],
            marker="o",
            linewidth=2,
            markersize=5,
            color=color,
        )

for ax in (ax1, ax2):
    ax.set_xticks(positions, labels)
    ax.set_xlim(-0.35, 5.85)
    ax.axvline(2.75, color="0.7", linewidth=0.8)
    ax.grid(axis="y", alpha=0.2)
    ax.set_axisbelow(True)
    ax.set_xlabel("Umap track and validity threshold")

ax1.set_title(
    "A  |  Matched narrowing score",
    fontsize=11,
    pad=12,
)
ax1.set_ylabel("log₂(reference W80 / peak W80)")
ax1.set_ylim(-0.025, 0.48)

ax2.set_title(
    "B  |  Eligible peak pairs excluded",
    fontsize=11,
    pad=12,
)
ax2.set_ylabel("Eligible pairs removed (%)")
ax2.set_ylim(0, 50)

fig.suptitle(
    "HCT116 narrowing sensitivity to mappability masking",
    fontsize=14,
    y=0.975,
)

handles, legend_labels = ax1.get_legend_handles_labels()
fig.legend(
    handles,
    legend_labels,
    loc="upper center",
    bbox_to_anchor=(0.5, 0.91),
    ncol=2,
    frameon=False,
)

fig.subplots_adjust(
    left=0.09,
    right=0.97,
    top=0.79,
    bottom=0.27,
    wspace=0.30,
)

fig.text(
    0.5, 0.105,
    "Dotted lines: unmasked scores | "
    "Pair loss measures eligible pairs, not sequencing reads",
    ha="center",
    fontsize=8.7,
)

fig.text(
    0.5, 0.065,
    "Umap k50/k100 sensitivity settings | "
    "Read length unknown | Exploratory; QC unassessed",
    ha="center",
    fontsize=8.7,
)

out = Path("figures/stage1")
out.mkdir(parents=True, exist_ok=True)

stem = out / "hct116_chr10_mappability_sensitivity"

fig.savefig(stem.with_suffix(".pdf"))
fig.savefig(stem.with_suffix(".png"), dpi=300)
plt.close(fig)

with Image.open(stem.with_suffix(".png")) as im:
    im.convert("RGB").save(
        out / "hct116_chr10_mappability_sensitivity_gallery.jpg",
        quality=95,
    )

print("PASS: Mappability result integrity verified")
print("PASS: All six masking scenarios plotted")
print("PASS: Mappability sensitivity PDF created")
print("PASS: High-resolution PNG created")
print("PASS: Gallery-compatible JPEG created")
