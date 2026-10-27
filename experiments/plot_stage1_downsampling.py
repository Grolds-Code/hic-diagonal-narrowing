"""Plot HCT116 W80 sensitivity to conditional contact thinning."""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image

base = Path("research_artifacts/stage1")

df = pd.read_csv(base / "w80_downsampling_v1.csv")
meta = json.loads(
    (base / "w80_downsampling_v1.metadata.json").read_text()
)

retentions = [0.50, 0.25, 0.10, 0.05, 0.01]
summaries = meta["summaries"]

assert len(df) == 610
assert len(summaries) == len(retentions)
assert not df.isna().any().any()

groups = []

for r, summary in zip(retentions, summaries):
    group = df[np.isclose(df.retention, r)].copy()

    assert len(group) == 122
    assert np.isclose(summary["retention"], r)
    assert group.replicates.eq(100).all()
    assert group.valid_replicates.eq(100).all()

    assert np.isclose(
        group.fraction_exact.mean(),
        summary["fraction_exact"],
    )

    groups.append(group)

print("PASS: All five retention levels verified")

x = np.arange(len(retentions))
labels = ["50%", "25%", "10%", "5%", "1%"]

exact = [
    100 * s["fraction_exact"] for s in summaries
]

within_one = [
    100 * s["fraction_within_one_bin"] for s in summaries
]

spans = [
    (
        g.thinned_w80_q95_bins
        - g.thinned_w80_q05_bins
    ).to_numpy() * 0.5
    for g in groups
]

fig, (ax1, ax2) = plt.subplots(
    1, 2,
    figsize=(10.5, 5.5),
)

ax1.plot(
    x, exact,
    color="#315D84",
    marker="o",
    linewidth=2,
    markersize=5,
    label="Exact W80 agreement",
)

ax1.plot(
    x, within_one,
    color="#438574",
    marker="s",
    linewidth=2,
    markersize=5,
    label="Within ±500 bp",
)

ax1.set_xticks(x, labels)
ax1.set_xlim(-0.25, 4.25)
ax1.set_ylim(0, 105)

ax1.set_xlabel("Contacts retained")
ax1.set_ylabel("Agreement with original W80 (%)")
ax1.set_title(
    "A  |  W80 agreement",
    fontsize=11,
    pad=12,
)

ax1.legend(
    loc="lower left",
    frameon=False,
    fontsize=9,
)
ax1.grid(axis="y", alpha=0.2)

bp = ax2.boxplot(
    spans,
    positions=x,
    widths=0.5,
    patch_artist=True,
    showfliers=True,
    medianprops={
        "color": "#263746",
        "linewidth": 2,
    },
    flierprops={
        "marker": ".",
        "markersize": 3,
        "alpha": 0.4,
    },
)

for patch in bp["boxes"]:
    patch.set_facecolor("#B8C6D6")

ax2.set_xticks(x, labels)
ax2.set_xlabel("Contacts retained")
ax2.set_ylabel("Per-window 90% thinning span (kb)")
ax2.set_title(
    "B  |  W80 spread across replicates",
    fontsize=11,
    pad=12,
)
ax2.grid(axis="y", alpha=0.2)

for ax in (ax1, ax2):
    ax.set_axisbelow(True)

fig.suptitle(
    "HCT116 W80 sensitivity to contact downsampling",
    fontsize=14,
    y=0.975,
)

fig.subplots_adjust(
    left=0.09,
    right=0.97,
    bottom=0.25,
    top=0.83,
    wspace=0.32,
)

fig.text(
    0.5, 0.085,
    "122 overlapping windows | 100 thinning replicates per depth",
    ha="center",
    fontsize=9,
)

fig.text(
    0.5, 0.045,
    "90% span = q95 − q05; not a confidence interval "
    "or biological replication",
    ha="center",
    fontsize=9,
)

out = Path("figures/stage1")
out.mkdir(parents=True, exist_ok=True)

stem = out / "hct116_chr10_w80_downsampling"

fig.savefig(stem.with_suffix(".pdf"))
fig.savefig(stem.with_suffix(".png"), dpi=300)
plt.close(fig)

with Image.open(stem.with_suffix(".png")) as im:
    im.convert("RGB").save(
        out / "hct116_chr10_w80_downsampling_gallery.jpg",
        quality=95,
    )

print("PASS: W80 downsampling PDF created")
print("PASS: High-resolution PNG created")
print("PASS: Gallery-compatible JPEG created")
