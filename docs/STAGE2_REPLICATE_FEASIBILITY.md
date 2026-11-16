# Stage 2 restart: ENCODE replicate feasibility

**Status:** Exploratory metadata audit, not a biological replication result.

## Purpose

Stage 1 identified five candidate narrowing footprints in HCT116 chr10:20–21 Mb
at 500-bp resolution. Coverage screening provisionally retained R1 and R2,
while R3–R5 were withheld because of spatial contact-support warnings.
The retained signals are sensitive to Umap masking, and no formal biological
validation has been established. See
[`STAGE1_SUPERVISOR_QC_UPDATE.md`](STAGE1_SUPERVISOR_QC_UPDATE.md).

We are resuming Stage 2 while awaiting Dr. Shamim's feedback, without
altering the published Stage 1 results.

## What this branch actually does

1. Query ENCODE experiment `ENCSR477GZK` and its released-file metadata.
2. Identify processed Hi-C contact matrices and BAM / pairs products.
3. Fetch detailed metadata for relevant processed products, including
   `derived_from` links and replicate-number fields when available.
4. Save a structured JSON inventory as a **GitHub Actions artifact**.
5. Test the metadata classification logic using offline unit tests.

No genomic data (FASTQ, BAM, Hi-C or contact counts) are downloaded. The
workflow is read-only with respect to ENCODE and does not push changes to
the repository.

Run the audit via the draft PR's GitHub Actions checks, or manually:

```bash
python -m unittest discover -s tests -p "test_stage2_encode_inventory.py" -v
python experiments/stage2_encode_replicate_inventory.py \
  --output research_artifacts/stage2/encode_replicate_inventory_v1.json
```

The command refuses to overwrite an existing output.

## Interpretation guardrails

- Replicate labels alone do not establish that a contact matrix is
  independently processed or suitable for biological-replicate validation.
- A matrix derived from multiple alignment files may be pooled.
- A replicate-specific BAM is not interchangeable with a processed,
  normalization-compatible replicate-specific contact matrix.
- Isogenic replicate groups may test reproducibility under shared genetic
  background, not biological generalization across independent donors or lines.
- R1 and R2 remain **provisional**; no selected candidate is confirmed as
  biologically significant.
- The earlier Stage 2 candidate/control aggregates were drawn from the same
  discovery matrix and are descriptive, not independent validation.

## Work intentionally not included in this PR

Local untracked artifacts from the previous Stage 2 exploration are **not**
present on GitHub: `modal_hct116_stage2_aggregate.py`,
`modal_hct116_stage2_umap_audit.py`, `stage2_prepare_aggregate_manifest.py`,
`stage2_plot_aggregates.py`, `figures/stage2/`, and
`research_artifacts/stage2/`. This branch does not reconstruct, overwrite,
or claim to validate those files.

## Next decision after metadata inspection

- **Independent replicate matrices available:** Check assay provenance,
  normalization, binning, per-replicate coverage, candidate reproducibility,
  and blind candidate/control comparison.
- **Only pooled matrices available:** Evaluate whether the linked alignment
  files permit separate replicate contact extraction without changing the
  experimental construct.
- **No suitable separated data available:** Continue a clearly labeled
  exploratory, coverage- and mappability-matched analysis, and identify
  independent data needed for biological validation.

Do not merge this preliminary feasibility work into a claim of Stage 2
biological replication without further verification.
