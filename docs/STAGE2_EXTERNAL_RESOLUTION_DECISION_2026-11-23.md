# Stage 2 external validation: resolution feasibility decision

**Decision recorded:** 2026-11-23. **Source:** user-run immutable local
`research_artifacts/stage2/hct116_external_hic_resolution_audit_v2.json`
(not checked into GitHub). **Status:** resolution metadata only; no local
contact blocks, W80, or biological concordance have been examined.

## Observed metadata

All four control HCT116 WT, untreated 0h sample matrices from GEO
GSE310051 were successfully parsed as HIC v8 by the bounded v2 audit.

| Sample | BP resolutions | chr10 name | chr10 length |
| --- | --- | --- | --- |
| GSM9287304 | 2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000, 1000 | chr10 | 133797422 |
| GSM9287305 | 2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000, 1000 | chr10 | 133797422 |
| GSM9287306 | 2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000, 1000 | chr10 | 133797422 |
| GSM9287307 | 2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000, 1000 | chr10 | 133797422 |

All four contain chr10:20,000,000–21,000,000 within the chromosome
extent. None advertises 500 bp in the BP resolution list. The finest
advertised BP resolution is 1000 bp.

**Assembly caveat:** The .hic genome field is a path
`/proj/phanstiel_lab/References/chromSizes/hg38_chromSizes.txt`.
This is consistent with hg38 labeling but did not pass the existing
exact-string genome check; verify against an authoritative chromosome-size
reference before using the maps for biological comparisons. The observed
chr10 length 133797422 is consistent with GRCh38.

## Scientific decision

1. **Do not** treat these four 1-kb maps as 500-bp replication evidence.
   Do not interpolate or upsample them to 500 bp.
2. Keep the original Stage 1 500-bp definition, 32-kb midpoint-assigned
   windows, maximum 32-kb separation, and predetermined regions frozen.
3. **Primary route:** locate external independent HCT116 replicates with
   actual 500-bp contact matrices and adequate local coverage, preferably
   with comparable experimental protocols.
4. **Secondary route:** consider an explicitly preregistered 1-kb
   sensitivity study on these four maps, separately labeled as
   cross-resolution and cross-protocol evidence. Confirm local raw
   contact coverage before proceeding; define valid-bin and depth/QC
   thresholds independently of W80 outcomes.
5. Preserve audit v1 and v2 locally, including error histories and
   timestamps. The header-only result is **not** a biological replication
   success or failure.

## Next action

Search source repositories and verified sample metadata for true
replicate-specific 500-bp HCT116 maps. Do not initiate full downloads or
any W80 biological analysis until source, independence, assembly,
resolution and coverage are verified.
