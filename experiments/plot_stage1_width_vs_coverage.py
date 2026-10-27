"""Plot aligned HCT116 W80 and contact-occupancy profiles."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pandas as pd
from PIL import Image

df = pd.read_csv(
    "research_artifacts/stage1/w80_profile_v2_offdiag.csv"
).sort_values("start_bp")

assert len(df) == 122
assert df.qc_status.eq("unassessed").all()

x = (df.start_bp + df.end_bp) / 2 / 1e6

fig, axes = plt.subplots(
    2, 1,
    figsize=(9, 7),
    sharex=True,
    gridspec_kw={"hspace": 0.15},
)

ax1, ax2 = axes

ax1.plot(x, df.w80_bp / 1000, color="#286092", lw=1.6, label="Local W80")
ax1.axhline(
    df.w80_bp.median() / 1000,
    color="#286092",
    ls="--",
    lw=1.1,
    alpha=0.7,
    label=f"Region-wide median ({df.w80_bp.median()/1000:.1f} kb)",
)
ax1.set_ylabel("W80 (kb)")
ax1.set_ylim(6, 10.5)
ax1.legend(loc="lower left", frameon=False, fontsize=9)

ax2.plot(
    x, df.pair_occupancy * 100,
    color="#386D5B",
    lw=1.7,
)
ax2.set_ylabel("Observed pair occupancy (%)")
ax2.set_xlabel("Genomic position on chr10 (Mb; GRCh38)")
ax2.set_ylim(0, 60)

candidates = [
    (20.296, 20.336, "A", "#8293A7"),
    (20.720, 20.792, "B", "#BA9770"),
]

for ax in axes:
    for start, end, name, color in candidates:
        ax.axvspan(start, end, color=color, alpha=0.20)
    ax.grid(axis="y", alpha=0.20)
    ax.set_xlim(20, 21)

fig.suptitle(
    "HCT116 diagonal width and contact occupancy",
    fontsize=15,
    y=0.975,
)

handles = [
    Patch(facecolor=color, alpha=0.35, label=f"Candidate {name}")
    for _, _, name, color in candidates
]

fig.legend(
    handles=handles,
    loc="lower center",
    bbox_to_anchor=(0.5, 0.105),
    ncol=2,
    frameon=False,
)

fig.subplots_adjust(
    left=0.12,
    right=0.97,
    top=0.88,
    bottom=0.24,
)

fig.text(
    0.5, 0.045,
    "500 bp resolution | 32 kb windows | 8 kb step | "
    "Exploratory; QC unassessed",
    ha="center",
    fontsize=9,
)

out = Path("figures/stage1")
out.mkdir(parents=True, exist_ok=True)
stem = out / "hct116_chr10_width_vs_coverage"

fig.savefig(stem.with_suffix(".pdf"))
fig.savefig(stem.with_suffix(".png"), dpi=300)
plt.close(fig)

with Image.open(stem.with_suffix(".png")) as im:
    im.convert("RGB").save(
        out / "hct116_chr10_width_vs_coverage_gallery.jpg",
        quality=95,
    )

print("PASS: Aligned W80 and occupancy PDF created")
print("PASS: High-resolution PNG created")
print("PASS: Gallery-compatible JPEG created")
