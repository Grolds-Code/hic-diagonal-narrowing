"""Plot descriptive HCT116 contact-pair occupancy."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
from PIL import Image

source = Path(
    "research_artifacts/stage1/w80_profile_v2_offdiag.csv"
)
output = Path("figures/stage1")
output.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(source).sort_values("start_bp")

assert len(df) == 122
assert df.qc_status.eq("unassessed").all()

x = (df.start_bp + df.end_bp) / 2 / 1e6
y = 100 * df.pair_occupancy

fig, ax = plt.subplots(figsize=(9, 4.5))

ax.plot(x, y, color="#294F78", linewidth=1.7)

ax.set_xlim(20, 21)
ax.set_ylim(0, 60)

ax.set_xlabel("Genomic position on chr10 (Mb; GRCh38)")
ax.set_ylabel("Observed contact-pair occupancy (%)")
ax.set_title(
    "HCT116 local Hi-C contact occupancy",
    fontsize=14,
    pad=15,
)

ax.grid(axis="y", alpha=0.2)

fig.subplots_adjust(
    left=0.11,
    right=0.97,
    bottom=0.22,
    top=0.86,
)

fig.text(
    0.5, 0.055,
    "500 bp resolution | 32 kb windows | 8 kb step | "
    "QC unassessed",
    ha="center",
    fontsize=9,
)

stem = output / "hct116_chr10_pair_occupancy"

fig.savefig(stem.with_suffix(".pdf"))
fig.savefig(stem.with_suffix(".png"), dpi=300)

plt.close(fig)

with Image.open(stem.with_suffix(".png")) as im:
    im.convert("RGB").save(
        stem.with_name(stem.name + "_gallery.jpg"),
        quality=95,
    )

print("PASS: Contact-occupancy PDF created")
print("PASS: Contact-occupancy PNG created")
print("PASS: Gallery-compatible JPEG created")
