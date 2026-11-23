# Stage 2: external HCT116 HIC v8 resolution metadata preflight

This audit is the **next step after** `stage2_external_hic_range_check.py`.
It reads only the first 65,536 bytes from each of four GEO HCT116 WT
untreated control .hic files, using exact HTTP 206 Content-Range responses.
It does not download complete files or read matrix blocks.

## Run

From the WSL research checkout, branch `research/stage2-replicate-audit`:

```bash
git pull --ff-only
uv run --with pytest python -m pytest -q
uv run python -m experiments.stage2_external_hic_resolution_audit
```

The output is written **once** to
`research_artifacts/stage2/hct116_external_hic_resolution_audit_v1.json`.
The script refuses to overwrite an existing artifact. Preserve that file and
its timestamp for provenance.

## Interpretation

A sample passes the *metadata gate* only if:
- the file parses as HIC version 8;
- the embedded reference identifier is `hg38` or `GRCh38`;
- chromosome 10 exists and contains chr10:20,000,000–21,000,000; and
- the BP zoom resolution list includes 500.

A failed range request or malformed header is reported as unverified, not
as evidence that the sample biologically lacks 500-bp data.

**Passing the gate does not demonstrate sufficient local depth.**
The script does **not** inspect any observed contact blocks, estimate
occupancy, missing-bin fraction, marginal read support, W80, or replication.
A local-coverage audit with a validated bounded contact extractor is the next
necessary gate before biological analysis. Do not tune the frozen Stage 1
estimator or move candidate loci after viewing validation results.

Maximum requested payload: 65,536 bytes per sample, 262,144 bytes in total.
The parser rejects unexpected HIC versions and invalid metadata counts.
If the server ignores HTTP Range, the script refuses to read its body.
