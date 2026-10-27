"""Plot the verified Stage 1 HCT116 W80 profile."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd

source = Path(
    "research_artifacts/stage1/w80_profile_v2_offdiag.csv"
)
output = Path("figures/stage1")
output.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(source).sort_values("start_bp")

assert len(df) == 122
assert df.qc_status.eq("unassessed").all()

x = (df.start_bp + df.end_bp) / 2 / 1_000_000
y = df.w80_bp / 1000
median = y.median()

fig, ax = plt.subplots(figsize=(9, 4.8))

ax.plot(x, y, linewidth=1.7, label="Local W80")
ax.axhline(
    median, linestyle="--", linewidth=1.2,
    label=f"Region-wide median = {median:.1f} kb"
)

ax.set(
    xlim=(20, 21),
    xlabel="Genomic position on chr10 (Mb; GRCh38)",
    ylabel="W80 (kb)",

)


# Exploratory candidate footprints, not biological boundaries.
candidates = [
    (20.296, 20.336, "A"),
    (20.720, 20.792, "B"),
]

for start, end, name in candidates:
    fill = "#8293A7" if name == "A" else "#BA9770"
    ax.axvspan(
        start, end,
        color=fill,
        alpha=0.18,
        label=f"Exploratory candidate {name}"
    )

ax.grid(alpha=0.2)
ax.legend(
    loc="lower center",
    bbox_to_anchor=(0.5, 1.025),
    ncol=2,
    frameon=False,
    fontsize=8.5,
)

fig.text(
    0.5, 0.025,
    "500 bp resolution | 32 kb windows | 8 kb step | "
    "Dmax = 32 kb | Exploratory; QC unassessed",
    ha="center", fontsize=8
)

fig.subplots_adjust(
    left=0.09,
    right=0.985,
    bottom=0.23,
    top=0.74,
)
fig.suptitle(
    "HCT116 local Hi-C diagonal-width profile",
    fontsize=15,
    y=0.975,
)

for extension in ("pdf", "png"):
    path = output / f"hct116_chr10_w80_profile.{extension}"
    fig.savefig(path, dpi=300)
    assert path.stat().st_size > 0
    print(f"PASS: {path}")

plt.close(fig)
