"""Plot conditional-thinning candidate-detection stability."""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image

base = Path("research_artifacts/stage1")
df = pd.read_csv(base / "candidate_stability_v1.csv")
meta = json.loads(
    (base / "candidate_stability_v1.metadata.json").read_text()
)

retentions = [0.25, 0.10]

assert len(df) == 200
assert not df.isna().any().any()
assert not df.duplicated(["retention", "replicate"]).any()
assert set(df.retention) == set(retentions)

groups = [
    df[np.isclose(df.retention, r)].sort_values("replicate")
    for r in retentions
]

for group, summary in zip(groups, meta["summaries"]):
    assert len(group) == meta["replicates_per_depth"] == 100
    assert np.isclose(group.retention.iloc[0], summary["retention"])
    assert np.isclose(
        group.region_a_retained.mean(),
        summary["region_persistence"][0],
    )
    assert np.isclose(
        group.region_b_retained.mean(),
        summary["region_persistence"][1],
    )
    assert np.isclose(
        group.newly_called_windows.median(),
        summary["median_newly_called_windows"],
    )

print("PASS: CSV results agree with experiment metadata")

a = [100 * g.region_a_retained.mean() for g in groups]
b = [100 * g.region_b_retained.mean() for g in groups]
new_calls = [g.newly_called_windows.to_numpy() for g in groups]

fig, (ax1, ax2) = plt.subplots(
    1, 2, figsize=(10, 5.4)
)

x = np.arange(len(retentions))
width = 0.32

bars_a = ax1.bar(
    x - width / 2, a, width,
    color="#48668A", label="Region A"
)
bars_b = ax1.bar(
    x + width / 2, b, width,
    color="#B17D48", label="Region B"
)

ax1.bar_label(bars_a, fmt="%.0f%%", padding=3, fontsize=9)
ax1.bar_label(bars_b, fmt="%.0f%%", padding=3, fontsize=9)

ax1.set_ylim(0, 115)
ax1.set_xticks(x, ["25%", "10%"])
ax1.set_xlabel("Contacts retained")
ax1.set_ylabel("Replicates retaining region (%)")
ax1.set_title("A  |  Candidate persistence", fontsize=11, pad=12)
ax1.grid(axis="y", alpha=0.2)
ax1.set_axisbelow(True)

ax1.legend(
    loc="upper center",
    bbox_to_anchor=(0.5, -0.20),
    ncol=2,
    frameon=False,
    fontsize=9,
)

bp = ax2.boxplot(
    new_calls,
    positions=[1, 2],
    widths=0.5,
    patch_artist=True,
    showfliers=True,
    flierprops={
        "marker": ".",
        "markersize": 3,
        "alpha": 0.45,
    },
    medianprops={"color": "#263746", "linewidth": 2},
)

for patch, color in zip(
    bp["boxes"], ["#AFC0D1", "#CDB79B"]
):
    patch.set_facecolor(color)

ax2.set_xticks([1, 2], ["25%", "10%"])
ax2.set_xlabel("Contacts retained")
ax2.set_ylabel("Newly called windows per replicate")
ax2.set_title("B  |  Additional detections", fontsize=11, pad=12)
ax2.grid(axis="y", alpha=0.2)
ax2.set_axisbelow(True)

fig.suptitle(
    "HCT116 candidate stability under contact thinning",
    fontsize=14,
    y=0.975,
)

fig.subplots_adjust(
    left=0.09,
    right=0.97,
    bottom=0.27,
    top=0.83,
    wspace=0.31,
)

fig.text(
    0.5, 0.055,
    "100 replicates per depth | Exploratory threshold = 0.2 | "
    "Conditional thinning; not biological replication",
    ha="center",
    fontsize=8.5,
)

out = Path("figures/stage1")
out.mkdir(parents=True, exist_ok=True)

stem = out / "hct116_chr10_candidate_stability"
fig.savefig(stem.with_suffix(".pdf"))
fig.savefig(stem.with_suffix(".png"), dpi=300)
plt.close(fig)

with Image.open(stem.with_suffix(".png")) as im:
    im.convert("RGB").save(
        out / "hct116_chr10_candidate_stability_gallery.jpg",
        quality=95,
    )

print("PASS: Candidate-stability PDF created")
print("PASS: High-resolution PNG created")
print("PASS: Gallery-compatible JPEG created")
print("NOTE: Newly called windows are not necessarily false positives.")
