# Stage 2 restart: Modal replicate-data and pooled-contact preflight

**Status:** Runnable Modal research workflow. Not a biological-replication
result. Dr. Shamim's Stage 1 QC feedback remains pending.

## Purpose

Stage 1 found five exploratory HCT116 narrowing footprints at 500-bp
resolution (GRCh38 chr10:20–21 Mb). A provisional coverage screen retains
R1 and R2, but mappability changes the calls. R3–R5 have low-support flags.
See [the Stage 1 supervisor update](STAGE1_SUPERVISOR_QC_UPDATE.md).

Stage 2 needs a credible candidate/control aggregate comparison and,
critically, independently processed biological-replicate evidence. We can
prepare that while retaining the existing untracked local Stage 2 work.

## Modal implementation

Run from a checkout with existing Modal credentials:

```bash
modal run modal_hct116_stage2_replicate_preflight.py
```

This uses Modal **CPU**, not a GPU. It attaches the existing
`hic-diagonal-narrowing-data` Volume at `/research` and does **not**
write to it. It reads:

- `/research/matrices/hct116_chr10_20000000_21000000_500bp_midpoint_32kb.npz`
- `/research/metadata/hct116_chr10_20000000_21000000_500bp_midpoint_32kb.json`
- The version-controlled Stage 1
  `research_artifacts/stage1/hct116_stage1_coverage_triage_v1.csv`

The remote function:

1. Verifies the extraction SHA256 against the pinned Stage 1 source and
   its provenance metadata.
2. Recomputes all 1,937 original midpoint-assigned W80 measurements
   (32-kb windows at 500-bp steps) and checks their contact totals against
   the archived Stage 1 table.
3. Confirms the 93 original candidate windows and 25 provisional
   coverage-screened calls (23 R1, 2 R2).
4. Queries ENCODE `ENCSR477GZK` public **metadata only** and checks
   processed contact products, individual BAM accessions, replicate
   labels and source lineage.
5. Returns JSON to the local Modal entrypoint, which saves it as
   `research_artifacts/stage2/hct116_stage2_modal_replicate_preflight_v1.json`.

The output filename is deliberately new. The script fails if it exists
already. Do not rename existing Stage 2 results or overwrite them.

**No BAM, FASTQ, .pairs or .hic file is downloaded.** Only the existing
small pooled 500-bp extracted .npz is read from Modal, plus remote ENCODE
metadata records. There is no expensive GPU allocation or background
deployment in this step.

## Tests

```bash
uv run --with numpy --with pytest pytest -q \
  tests/test_stage2_encode_inventory.py \
  tests/test_stage2_pooled_preflight.py
```

Modal execution is not authenticated through the GitHub connection.
An authenticated Modal CLI is required to start the run. The JSON
report is created only after a successful remote computation. There is
no GitHub Actions workflow on this branch.

## Findings already established

A previous **metadata-only GitHub audit** (before switching to Modal)
identified four .hic records and two separate biological-replicate BAMs.
The two matrices with explicit lineage each derive from **both** BAMs;
their pooled-source W80 cannot establish biological replication.
Other .hic products lack sufficient lineage metadata. Historical
provenance remains in
[STAGE2_REPLICATE_FEASIBILITY_FINDINGS.md](STAGE2_REPLICATE_FEASIBILITY_FINDINGS.md).
The historical GitHub Actions run is not the execution route for
future analyses.

## Guardrails

- A copied 500-bp W80 result verifies only software and source provenance.
- R1 and R2 remain exploratory, not validated.
- Overlapping 32-kb windows are not independent measurements.
- Two biological-replicate BAMs do not equal two independent Hi-C maps.
- Biological replicate and control comparisons require
  comparable processing, coverage and mappability checks.
- The initial Stage 2 aggregates are descriptive and reside in local
  **untracked** files. They were not added, modified or reproduced in
  this PR.
- Do not call a copy-number deletion or a biological mechanism
  based solely on contact dropouts.

## Next computational milestone

First establish the schema and data-provenance requirements for deriving
**two separately processed 500-bp contact maps** from
`ENCFF624NZV` and `ENCFF583BZD`, without converting mapped Hi-C
alignments into contact pairs using unjustified assumptions.

Then, only after comparable maps and pre-specified controls are
available, evaluate R1/R2 W80 reproducibility and error modes per
replicate. The 200-bp and 100-bp analyses remain separate, pending
coverage/resolution feasibility checks.
