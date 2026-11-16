# Stage 2: bounded replicate-BAM schema probe on Modal

**Status:** Implemented on the Stage 2 research branch; execution and
unit-test outcomes still pending.

**Motivation:** The previously verified ENCODE metadata audit identified
biological-replicate BAMs `ENCFF624NZV` (replicate 1) and `ENCFF583BZD`
(replicate 2), but did not establish two independent .hic matrices.
The [pooled-data Modal preflight](STAGE2_MODAL_PREFLIGHT_RESULT.md)
successfully reproduced all 1,937 archived W80 measurements. That is a
software/provenance check, not replicate concordance.

## Low-cost probe

```bash
uv run --with pytest pytest -q \
  tests/test_stage2_encode_inventory.py \
  tests/test_stage2_pooled_preflight.py \
  tests/test_stage2_bam_bounded_schema.py

uv run --with modal modal run \
  modal_hct116_stage2_bam_schema_probe.py
```

The Modal app uses **CPU only** and reads ENCODE metadata, then attempts
a strict HTTP byte-range request for the **first at most 2 MiB of
each BAM** (4 MiB across two files). It refuses to stream the file when
the server does not respond with HTTP 206 and a range starting at byte
zero. It decodes only *complete* BGZF blocks from the retrieved prefix,
and reports BAM header sequence dictionary / sort-order, and up to
2,000 *initial* alignment record flags and MAPQ values.

The output is local JSON:
`research_artifacts/stage2/hct116_stage2_bam_bounded_schema_v1.json`.
Existing files are never overwritten; the Modal job does not mount or
modify the existing HCT116 volume.

## Interpreting results

- `bounded_prefix_obtained` means a small BAM prefix was obtained and
  parsed. **It does not prove that valid Hi-C ligation-contact pairs
  can be reconstructed** from BAM SAM mate fields.
- A single initial coordinate-sorted BAM prefix is **not** representative
  of the rest of the file. Its flag and MAPQ counts must not be used
  for biological interpretation, inferential QC or read-depth estimates.
- `range_not_supported` or a bounded-probe error means the format
  sample was inconclusive; it is **not** evidence that the replicate
  contains no usable contact pairs.
- The data source is ENCODE's public release; the read headers may
  identify pipeline sorting / reference dictionaries, not Hi-C
  contact-pair semantics or duplicate policy.
- The 500-bp aggregate Stage 2 comparison and true biological-replicate
  concordance are **still outstanding**.
- We deliberately do **not** download whole BAMs or allocate GPUs at
  this stage.

## Decision gate after the run

If both headers can be parsed, inspect the alignment flags and ENCODE
pipeline documentation before designing any conversion from BAM
alignments into independent contact pairs. If the semantics are
unclear, do **not** fabricate contact matrices from BAM mate positions.
Look for documented ENCODE post-processing steps, replicate-specific
contacts or processing output instead.

Next, derive or obtain two genuinely comparable, separate
500-bp matrices and test coverage, mappability, W80 and fixed
candidate/control comparisons per replicate. Keep R1/R2 as
provisional until that validation is complete.
