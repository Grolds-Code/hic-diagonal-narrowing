"""Visual comparison of Stage 1 coverage screening and Umap sensitivity."""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

root = Path("research_artifacts/stage1")
figure_path = Path(
    "figures/stage1/"
    "hct116_stage1_supervisor_qc_comparison_v1.png"
)

if figure_path.exists():
    raise RuntimeError(
        f"Refusing to overwrite: {figure_path}"
    )

regions = pd.read_csv(
    root / "hct116_stage1_candidate_deliverable_v1.csv"
)
triage = pd.read_csv(
    root / "hct116_stage1_coverage_triage_v1.csv"
)
umap = pd.read_csv(
    root / "hct116_stage1_umap_full_rescore_v1.csv"
)

if (
    len(regions) != 5
    or len(triage) != 1937
    or len(umap) != 1809 * 7
):
    raise RuntimeError("Unexpected Stage 1 dimensions")

def calls(series):
    return series.astype(str).str.lower().eq("true")

labels = []
original_counts = []
screened_counts = []

for region in regions.itertuples(index=False):
    labels.append(region.candidate_id)

    membership = (
        (triage["start_bp"] >= region.footprint_start_bp)
        & (triage["end_bp"] <= region.footprint_end_bp)
    )

    original_counts.append(int(
        (membership & calls(
            triage["original_candidate_call"]
        )).sum()
    ))

    screened_counts.append(int(
        (membership & calls(
            triage["preview_candidate_call"]
        )).sum()
    ))

if original_counts != [23, 2, 8, 11, 49]:
    raise RuntimeError("Original call mismatch")

if screened_counts != [23, 2, 0, 0, 0]:
    raise RuntimeError("Screened call mismatch")

scenarios = [
    "unmasked",
    "k50_50",
    "k50_80",
    "k50_95",
    "k100_50",
    "k100_80",
    "k100_95",
]

r1_counts = []
r2_counts = []
other_counts = []

for scenario in scenarios:
    subset = umap.loc[
        (umap["scenario"] == scenario)
        & calls(umap["candidate_call"])
    ]

    region_counts = {}

    for label in ("R1", "R2"):
        region = regions.loc[
            regions["candidate_id"] == label
        ].iloc[0]

        region_counts[label] = int((
            (subset["start_bp"] >= region.footprint_start_bp)
            & (subset["end_bp"] <= region.footprint_end_bp)
        ).sum())

    r1_counts.append(region_counts["R1"])
    r2_counts.append(region_counts["R2"])
    other_counts.append(
        len(subset)
        - region_counts["R1"]
        - region_counts["R2"]
    )

fig, axes = plt.subplots(
    2, 1,
    figsize=(11, 8),
    constrained_layout=True,
)

# Panel A: coverage-screening effect.
x = np.arange(len(labels))
width = 0.34

a = axes[0].bar(
    x - width / 2,
    original_counts,
    width,
    label="Original detector",
)

b = axes[0].bar(
    x + width / 2,
    screened_counts,
    width,
    label="Coverage-screened preview",
)

axes[0].bar_label(a, padding=3)
axes[0].bar_label(b, padding=3)
axes[0].set_xticks(x, labels)
axes[0].set_ylabel("Called 500-bp scan windows")
axes[0].set_title(
    "A. Original vs provisional coverage-screened calls"
)
axes[0].legend(frameon=False)
axes[0].set_ylim(0, 56)

# Panel B: mappability sensitivity.
x2 = np.arange(len(scenarios))

for name, counts in (
    ("R1 footprint", r1_counts),
    ("R2 footprint", r2_counts),
    ("Other positions", other_counts),
):
    axes[1].plot(
        x2,
        counts,
        marker="o",
        linewidth=1.8,
        label=name,
    )

axes[1].set_xticks(
    x2,
    [
        "Unmasked",
        "k50\n0.50",
        "k50\n0.80",
        "k50\n0.95",
        "k100\n0.50",
        "k100\n0.80",
        "k100\n0.95",
    ],
)

axes[1].set_ylabel("Called 500-bp scan windows")
axes[1].set_title(
    "B. Umap sensitivity after coverage screening"
)
axes[1].legend(frameon=False)
axes[1].grid(axis="y", alpha=0.2)

fig.suptitle(
    "HCT116 chr10:20–21 Mb | Stage 1 QC sensitivity",
    fontsize=14,
)

fig.text(
    0.5, 0.005,
    "Exploratory results. Overlapping windows are not independent. "
    "QC thresholds and read-length assumptions are provisional.",
    ha="center",
    fontsize=8,
)

figure_path.parent.mkdir(
    parents=True,
    exist_ok=True,
)

fig.savefig(
    figure_path,
    dpi=200,
    bbox_inches="tight",
)

plt.close(fig)

if not figure_path.is_file():
    raise RuntimeError("Figure not created")

print("\nSTAGE 1 — SUPERVISOR QC FIGURE")
print("=" * 65)
print("Original calls:", sum(original_counts))
print("Coverage-screened calls:", sum(screened_counts))
print("Umap scenarios:", len(scenarios))
print("R1 call counts:", r1_counts)
print("R2 call counts:", r2_counts)
print("Other call counts:", other_counts)
print("Saved:", figure_path)
print("Step 392: PASS")
