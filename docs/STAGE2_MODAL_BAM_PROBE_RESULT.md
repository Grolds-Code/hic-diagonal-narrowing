# Stage 2 Modal biological-replicate BAM schema probe: verified result

**Date:** 2026-11-17 (Africa/Nairobi)  
**Status:** BAM prefix inspection succeeded; not contact-pair reconstruction or biological validation.  
**Modal execution:** [ap-CGOQYlE3SkfiewWdL2cuJ4](https://modal.com/apps/groldotieno97/main/ap-CGOQYlE3SkfiewWdL2cuJ4)  
**Script:** `modal_hct116_stage2_bam_schema_probe.py`  
**Local output (not committed):** `research_artifacts/stage2/hct116_stage2_bam_bounded_schema_v1.json` in `~/projects/hic-stage2-modal`

## Exact console observations

| Read-only measurement on **initial** 2,000 BAM alignment records | Biological replicate 1 | Biological replicate 2 |
|---|---:|---:|
| ENCODE BAM accession | `ENCFF624NZV` | `ENCFF583BZD` |
| Download/probe status | `bounded_prefix_obtained` | `bounded_prefix_obtained` |
| BAM reference sequences in header | 195 | 195 |
| Complete alignment records inspected | 2,000 | 2,000 |
| SAM paired-read flag | 0 | 0 |
| SAM proper-pair flag | 0 | 0 |
| SAM unmapped flag | 444 | 488 |
| SAM mate-unmapped flag | 0 | 0 |
| SAM secondary-alignment flag | 1,042 | 1,011 |
| SAM duplicate flag | 0 | 0 |
| SAM supplementary-alignment flag | 0 | 0 |
| Alignment records mapped to chr10 | 86 | 48 |

The Modal log explicitly reports:
- `Full BAM files downloaded: NO`
- `Hi-C replicate W80 analysis: NOT PERFORMED`
- The Modal app completed normally.

**Provenance:** These observations are transcribed from the researcher's
completed Modal console log. The locally saved JSON has not been
separately inspected or committed to this repository.

## Meaning and limitation

Both BAMs are valid enough for the bounded BGZF/BAM decoder to parse
the headers and 2,000 complete alignment records. Neither **initial**
alignment sample contains a record with the conventional SAM
paired-read flag set.

This is evidence **against assuming** that conventional BAM mate
coordinates can directly recover Hi-C ligation pairs, but it does **not**
prove all records are unpaired or prove how this specific pipeline
represents and joins contacts. An initial prefix of a potentially
coordinate-sorted BAM is not a representative sample: do not
extrapolate unmapped, secondary, MAPQ or chr10 proportions to the
full alignment product.

ENCODE's public Hi-C documentation distinguishes processed BAM
alignments and the `.pairs` contact list. Its Hi-C pipeline involves
contact selection, filtering and deduplication before matrix
generation. Independent replicate contact matrices must preserve
correct read-pair assembly, mapping-quality rules, duplicate handling,
ligation-aware filtering and depth/provenance consistency.

References:
- https://www.encodeproject.org/hic/
- https://pmc.ncbi.nlm.nih.gov/articles/PMC10104020/
- https://github.com/ENCODE-DCC/hic-pipeline

## Next research decision

Inspect the full *local JSON summary* for BAM header sort order,
read-group/program header counts, and file metadata; then match these
against the ENCODE Juicer pipeline's exact alignment-to-contact
conversion steps.

Do **not** convert independently mapped alignments into contacts
by blindly reading `RNEXT`/`PNEXT` or changing the W80 estimator.
Before processing full replicate BAMs, identify an ENCODE-supported
replicate-specific contact pipeline, estimate input/compute costs and
reproduce the pooled contact extraction as a provenance check.

R1 and R2 remain provisional and mappability-sensitive. Stage 2
biological replication remains **NOT TESTED**.
