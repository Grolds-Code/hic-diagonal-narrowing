"""Compare HCT116 model-conditional intervals with technical thinning.

The full-depth 95% model interval and the 5th-95th percentile
range of thinned estimates have different interpretations.
Their overlap is descriptive, NOT a confidence coverage test.
"""

import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from statistics import median

ROOT = Path("research_artifacts/stage1")

THINNING = ROOT / "w80_downsampling_v1.csv"
THINNING_META = ROOT / "w80_downsampling_v1.metadata.json"

INTERVALS = ROOT / "hct116_conditional_w80_intervals_v1.csv"

STABILITY_META = ROOT / "candidate_stability_v1.metadata.json"
STABILITY_CSV = ROOT / "candidate_stability_v1.csv"

PROFILE_META = ROOT / "w80_profile_v2_offdiag.metadata.json"


def load_csv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def check_hash(path, expected):
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    if actual != expected:
        raise RuntimeError(f"Checksum mismatch: {path}")


thin_meta = json.loads(THINNING_META.read_text())
stability_meta = json.loads(STABILITY_META.read_text())
profile_meta = json.loads(PROFILE_META.read_text())

check_hash(THINNING, thin_meta["output_sha256"])
check_hash(STABILITY_CSV, stability_meta["output_sha256"])

for meta in (thin_meta, stability_meta):
    if meta["source_sha256"] != profile_meta["source_sha256"]:
        raise RuntimeError("Inconsistent source provenance")

if thin_meta["reference_sha256"] != profile_meta["output_sha256"]:
    raise RuntimeError("Inconsistent W80 reference")

if stability_meta["profile_sha256"] != profile_meta["output_sha256"]:
    raise RuntimeError("Inconsistent candidate reference")

interval_rows = load_csv(INTERVALS)
thinning_rows = load_csv(THINNING)

intervals = {}

for row in interval_rows:
    key = (int(row["start_bp"]), int(row["end_bp"]))

    if key in intervals:
        raise RuntimeError("Duplicate interval window")

    lo = int(row["lower_95_bins"])
    w80 = int(row["w80_bins"])
    hi = int(row["upper_95_bins"])

    if not 1 <= lo <= w80 <= hi <= 64:
        raise RuntimeError("Invalid model interval")

    intervals[key] = row

if len(intervals) != 122:
    raise RuntimeError("Expected 122 model intervals")

by_depth = defaultdict(list)
seen = set()

for row in thinning_rows:
    key = (int(row["start_bp"]), int(row["end_bp"]))
    retention = float(row["retention"])

    if key not in intervals:
        raise RuntimeError("Unknown thinned window")

    if (retention, key) in seen:
        raise RuntimeError("Duplicate thinning window")

    seen.add((retention, key))

    reference = intervals[key]

    if int(row["baseline_w80_bins"]) != int(
        reference["w80_bins"]
    ):
        raise RuntimeError("Baseline W80 disagreement")

    if float(row["baseline_contact_signal"]) != float(
        reference["total_contacts"]
    ):
        raise RuntimeError("Baseline contact-count disagreement")

    if (
        int(row["replicates"]) != 100
        or int(row["valid_replicates"]) != 100
    ):
        raise RuntimeError("Unexpected replicate count")

    q05 = float(row["thinned_w80_q05_bins"])
    q95 = float(row["thinned_w80_q95_bins"])
    thinned_median = float(row["median_thinned_w80_bins"])

    if not 1 <= q05 <= thinned_median <= q95 <= 64:
        raise RuntimeError("Invalid thinning quantiles")

    lo = int(reference["lower_95_bins"])
    hi = int(reference["upper_95_bins"])

    by_depth[retention].append({
        "key": key,
        "original_w80": int(reference["w80_bins"]),
        "model_lower": lo,
        "model_upper": hi,
        "model_span": hi - lo,
        "thinned_q05": q05,
        "thinned_q95": q95,
        "thinned_span": q95 - q05,
        "thinned_median": thinned_median,
        "thinning_envelope_contained": (
            lo <= q05 and q95 <= hi
        ),
        "thinning_median_inside": (
            lo <= thinned_median <= hi
        ),
    })

expected_depths = set(
    float(x) for x in thin_meta["retention_fractions"]
)

if set(by_depth) != expected_depths:
    raise RuntimeError("Unexpected retention levels")

summary = {}

print("\nHCT116 CONDITIONAL INTERVALS VS TECHNICAL THINNING")
print("=" * 69)

print(
    "Retention  Median model span  Median thinning 90% span"
    "  Envelope contained"
)

for retention in sorted(by_depth, reverse=True):
    rows = by_depth[retention]

    if len(rows) != 122:
        raise RuntimeError("Incomplete thinning windows")

    model_span = median(r["model_span"] for r in rows)
    thinning_span = median(r["thinned_span"] for r in rows)

    envelope_rate = sum(
        r["thinning_envelope_contained"] for r in rows
    ) / len(rows)

    median_inside_rate = sum(
        r["thinning_median_inside"] for r in rows
    ) / len(rows)

    print(
        f"{retention:>7.0%}"
        f"{model_span:>19.1f}"
        f"{thinning_span:>26.2f}"
        f"{envelope_rate:>20.1%}"
    )

    summary[str(retention)] = {
        "windows": len(rows),
        "median_full_depth_model_span_bins": model_span,
        "median_thinning_90pct_span_bins": thinning_span,
        "fraction_thinning_envelopes_inside_full_interval": (
            envelope_rate
        ),
        "fraction_thinning_medians_inside_full_interval": (
            median_inside_rate
        ),
    }

candidates = {
    "Region A": (20_296_000, 20_328_000),
    "Region B": (20_736_000, 20_768_000),
}

print("\nCANDIDATE WINDOW DETAILS")
print("-" * 69)

for name, key in candidates.items():
    ref = intervals[key]

    print(
        f"\n{name}: original W80={ref['w80_bins']}, "
        f"full-depth model interval="
        f"[{ref['lower_95_bins']}, {ref['upper_95_bins']}]"
    )

    for retention in sorted(by_depth, reverse=True):
        matches = [
            row for row in by_depth[retention]
            if row["key"] == key
        ]

        if len(matches) != 1:
            raise RuntimeError("Missing candidate thinning result")

        row = matches[0]

        print(
            f"  {retention:>3.0%} retained: "
            f"thinning q05-q95="
            f"[{row['thinned_q05']:g}, "
            f"{row['thinned_q95']:g}], "
            f"median={row['thinned_median']:g}"
        )

print("\nEXISTING CANDIDATE-CALL PERSISTENCE")

for item in stability_meta["summaries"]:
    retention = float(item["retention"])

    for case, rate in zip(
        ("Region A", "Region B"),
        item["region_persistence"],
    ):
        print(
            f"  {case}, retention {retention:.0%}: "
            f"{rate:.0%} region persistence"
        )

output = ROOT / "hct116_interval_thinning_comparison_v1.json"

result = {
    "analysis": "descriptive_full_depth_interval_vs_thinning",
    "thinning_source_sha256": thin_meta["source_sha256"],
    "comparison_by_retention": summary,
    "limitation": (
        "Thinning quantile envelopes condition on observed source "
        "counts; full-depth model confidence intervals concern "
        "a model-defined target. Containment is not CI coverage."
    ),
}

output.write_text(json.dumps(result, indent=2) + "\n")

print("\nSaved:", output)
print("Diagnostic only; no empirical confidence-coverage claim.")
