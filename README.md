# Hi-C Diagonal Narrowing

## A Framework for Local Contact-Width Analysis in Hi-C Data

**Version:** 0.1.0
**Status:** Exploratory computational research

This project develops and evaluates methods for detecting localized narrowing in Hi-C contact-distance profiles, with particular attention to sequencing coverage, genomic mappability, and technical artifacts.

No biological significance is assumed from detection alone.

## 1. Scientific Motivation

Hi-C contact matrices describe interactions between genomic loci. Contact frequency generally decreases with genomic separation, producing a prominent band near the matrix diagonal.

Local changes in diagonal width may reflect variation in chromatin organization. However, similar patterns can arise from sparse observations, sequencing coverage, low sequence mappability, or mapping-quality filtering.

This research investigates whether local diagonal narrowing can be measured reproducibly while identifying technical conditions that undermine its interpretation.

## 2. Methodological Framework

The proposed method estimates local diagonal width from the genomic-distance distribution of observed cis-contacts.

For each genomic window, contact pairs are assigned by their genomic midpoint. Main-diagonal contacts are excluded, and analysis is restricted to a specified maximum separation.

Mean contact intensity is calculated at each genomic separation using the number of valid eligible bin pairs, including those with zero observed contacts.

The primary width statistic, W80, is the smallest genomic separation at which the cumulative mean contact intensity reaches 80% of the total within the analyzed separation range.

**Current implementation:** The core `estimate_w80()` function assumes
a complete, symmetric, nonnegative contact matrix. It does not currently
accept missing-bin or mappability masks. Mask-aware analyses are
implemented separately in experimental scripts.

The function accepts half-open genomic windows expressed in bin indices
and returns W80 in bins. Multiplication by the bin resolution converts
the result into base pairs.

### 2.1 Mathematical Definitions

For a zero-based, half-open genomic window $W=[a,b)$ expressed in bin coordinates, eligible contact pairs at separation $d$ are defined as:

$$
P_W(d)=\left\{(i,j)\;\middle|\;0\le i<j<n,\quad j-i=d,\quad a\le\frac{i+j}{2}<b\right\}
$$

**where:**

- $n$ is the total number of genomic bins in the contact matrix.
- $i$ and $j$ are genomic bin indices, with $i<j$ to count each undirected pair once.
- $d=j-i$ is the genomic separation in bins.
- $a$ and $b$ define the half-open analysis window.
- $(i+j)/2$ is the contact pair's genomic midpoint.

The current unmasked implementation treats all such pairs as eligible. Mask-aware experiments additionally exclude pairs whose endpoints fail the specified validity mask.

Mean contact intensity:

$$
m_W(d) = \frac{\sum_{(i,j)\in P_W(d)} C(i,j)}{|P_W(d)|}
$$

**where:**

- $m_W(d)$ is the mean observed contact intensity at separation $d$ within genomic window $W$.
- $W$ is the genomic analysis window.
- $d$ is genomic separation measured in bins.
- $C(i,j)$ is the observed Hi-C contact count between bins $i$ and $j$.
- $P_W(d)$ is the set of eligible contact pairs at separation $d$, assigned to $W$ by their genomic midpoints.
- $|P_W(d)|$ is the number of eligible pairs, including pairs with zero observed contacts.

**Interpretation:** Dividing by the number of eligible pairs accounts for differences in pair availability across genomic separations.

Local diagonal width:

$$
F_W(D) = \frac{\sum_{d=1}^{D}m_W(d)}{\sum_{d=1}^{D_{\max}}m_W(d)}
$$

$$
W_{80}(W) = \min\left\{D : F_W(D) \geq 0.80\right\}
$$

**where:**

- $F_W(D)$ is the cumulative fraction of mean contact intensity up to separation $D$.
- $D$ is a candidate genomic separation, measured in bins.
- $D_{\max}$ is the maximum analyzed separation.
- $m_W(d)$ is the mean contact intensity at separation $d$.
- $W_{80}(W)$ is the smallest separation at which at least 80% of the total mean contact intensity has accumulated.

**Interpretation:** A smaller W80 indicates that observed contact intensity is more concentrated at shorter genomic separations.

**Assumptions:** The total mean contact intensity must be positive, and each analyzed separation must have eligible contact pairs. The estimator returns W80 in bins; multiplying by genomic resolution converts it to base pairs.

Reference-width estimation:

$$
W_{80}^{\mathrm{ref}}(W)
=
\operatorname{median}_{V\in\mathcal{R}(W)}
W_{80}(V)
$$

**where:**

- $\mathcal{R}(W)$ is the set of eligible reference windows.
- $V$ denotes a reference window on the same chromosome as $W$.
- $W_{80}(V)$ is the estimated width of reference window $V$.
- $W_{80}^{\mathrm{ref}}(W)$ is the median reference width.

In the current detector, reference windows must not overlap
the target window. They must have finite, positive W80
values, positive contact signal, and `qc_status = "unassessed"`.

At least three eligible reference windows are required by
default. Reference windows may overlap one another and are
not automatically matched for sequencing coverage.

Relative narrowing score:

$$
S(W) = \log_2\left(
\frac{W_{80}^{\mathrm{ref}}(W)}{W_{80}(W)}
\right)
$$

**where:**

- $S(W)$ is the relative narrowing score for genomic window $W$.
- $W_{80}(W)$ is the local diagonal width.
- $W_{80}^{\mathrm{ref}}(W)$ is the reference width.
- $\log_2$ denotes the base-2 logarithm.

**Interpretation:**

- $S(W) > 0$: local narrowing relative to the reference.
- $S(W) = 0$: equal local and reference widths.
- $S(W) < 0$: local widening relative to the reference.

**Assumptions:** Both widths must be positive and expressed
in the same units. The score is descriptive and does not
constitute a statistical significance measure.

## 3. Software Architecture

The research implementation is organized as a Python package under `src/hic_diagonal_narrowing/`.

| Module | Scientific responsibility |
| --- | --- |
| `width.py` | Estimation of local diagonal width (W80) |
| `scan.py` | Sliding-window analysis of contact-width profiles |
| `detect.py` | Identification of relative narrowing candidates |
| `qc.py` | Contact-support and quality-control diagnostics |
| `uncertainty.py` | Bootstrap-based uncertainty estimation |
| `synthetic.py` | Generation of synthetic Hi-C contact matrices |

The `experiments/` directory contains numerical validation and calibration studies.

Scripts prefixed with `modal_` support remote HCT116 data extraction, contact reconstruction, and technical sensitivity analyses.

The software is currently a research prototype. A stable public API, unified command-line interface, and comprehensive automated test suite remain under development.

### Sliding-Window Implementation

The `scan_local_width()` function uses zero-based, half-open
genomic coordinates. Its default configuration uses 64-bin
windows, 16-bin steps, and a maximum separation of 64 bins.

At 500 bp resolution, these correspond to 32 kb windows
and 8 kb scanning steps.

The scanner reports W80, total contact signal, pair occupancy,
and distance-support diagnostics.

Windows without observed contact signal receive a missing
W80 value and `qc_status = "no_signal"`.

Other informative windows receive `qc_status = "unassessed"`.

The scanner does not yet apply calibrated QC acceptance thresholds.
Candidate calling is handled separately by the detection module.

### Candidate-Detection Implementation

The `detect_narrowing()` function calculates descriptive scores
using the median W80 of informative windows on the same chromosome
that do not overlap the target window.

These reference windows are not necessarily mutually disjoint or
matched for contact coverage. Coverage-matched comparisons were
investigated separately in experimental analyses.

The function's default score threshold is 0.5. The exploratory
HCT116 candidate analysis used a threshold of 0.20.

Overlapping or directly adjacent candidate windows are merged
into detection footprints. Their boundaries are not inferred
biological boundaries.

The detector currently accepts windows with
`qc_status = "unassessed"` and does not enforce calibrated QC
acceptance thresholds.

Scores and candidate calls are descriptive, not statistical
significance estimates.

### Contact-Support Quality Control

The `summarize_local_support()` function reports:

- Total observed contact signal across eligible pairs.
- Number of eligible and nonzero contact pairs.
- Pair occupancy: nonzero pairs divided by eligible pairs.
- Eligible and supported genomic-distance bands.
- Distance-support fraction: supported bands divided by eligible bands.

Only unique upper-triangle, off-diagonal pairs are included,
using midpoint-based window assignment.

The current implementation assumes a complete, symmetric,
nonnegative contact matrix without explicit missing-bin masks.

These metrics characterize observed data support. They are not
calibrated biological confidence measures and do not independently
establish whether a narrowing candidate is reliable.

### Bootstrap Uncertainty Estimation

The `bootstrap_w80()` function estimates conditional sampling
variability using multinomial resampling of observed contact
counts across genomic-distance bands.

The default configuration uses 300 bootstrap realizations,
a 95% confidence level, and random seed 42.

Each realization preserves the total observed contact count.
W80 is recalculated using mean intensity per eligible pair
at each genomic separation.

The reported percentile intervals describe uncertainty
under this conditional sampling model.

They do not account for spatial contact dependencies,
missing-bin effects, mapping bias, or biological-replicate
variability. They are not comprehensive biological
confidence intervals.

### Installation and Verification

Requires Python 3.12+ and `uv`. Run `uv sync` from the repository root.

Verify the estimator using `uv run python -c "from hic_diagonal_narrowing.width import estimate_w80; print(estimate_w80.__name__)"`.

## 4. Synthetic Validation

The W80 estimator was evaluated using controlled synthetic Hi-C contact matrices representing three conditions:

1. Baseline distance-dependent contact decay.
2. Engineered narrowing of the contact diagonal.
3. Uniform reduction in contact coverage without changing the underlying contact-distance profile.

### 4.1 Numerical Results

Across 100 independent simulated realizations per condition:

| Condition | Mean W80 (bins) | Standard deviation |
| --- | ---: | ---: |
| Baseline | 27.64 | 0.48 |
| Engineered narrowing | 8.96 | 0.20 |
| Coverage reduction | 27.58 | 0.50 |

The estimator distinguished engineered narrowing from uniform coverage reduction under these synthetic conditions.

### 4.2 Robustness and Uncertainty

Additional experiments investigated sensitivity to contact depth, bootstrap uncertainty, and the reliability of width estimates under increasingly sparse observations.

At extremely low simulated coverage, false narrowing detections emerged, demonstrating the need for explicit quality-control and reliability criteria.

These experiments establish preliminary numerical validity under controlled assumptions. They do not establish general reliability across all realistic Hi-C technical artifacts or biological datasets.

### 4.3 Synthetic Data-Generating Mechanism

The `simulate_contact_map()` function generates symmetric,
integer-valued contact matrices using Poisson sampling.

The default configuration uses 256 genomic bins and a
target region spanning bins 96–160.

Baseline expected contact intensity declines exponentially
with genomic separation, using a decay scale of 18 bins.

The narrowing scenario reduces the regional decay scale
to 5 bins and rescales expected contact counts to preserve
the regional expected total across simulated off-diagonal
pairs. This preservation does not necessarily hold within
every truncated analysis distance range.

The coverage-reduction scenario multiplies regional expected
contact counts by 0.25 without changing the decay shape.

This is a uniform expected-intensity reduction, not an
explicit model of missing bins, mapping failures, or
spatially heterogeneous sequencing coverage.

These simulations provide controlled numerical test cases
rather than realistic biological models of chromatin folding.

## 5. Real-Data Analysis: HCT116 Hi-C

### 5.1 Dataset Provenance

The initial real-data investigation uses publicly available intact Hi-C data from the human colorectal carcinoma cell line HCT116.

| Property | Specification |
| --- | --- |
| Data repository | ENCODE |
| Experiment accession | ENCSR477GZK |
| Primary file accession | ENCFF573OPJ |
| Genome assembly | GRCh38 (hg38) |
| Assay | Intact Hi-C |
| Analyzed chromosome | chr10 |
| Genomic interval | 20,000,000–21,000,000 bp |
| Resolution | 500 bp |
| Maximum contact separation | 32 kb |
| Analysis window | 32 kb |
| Window step | 8 kb |

The primary contact matrix combines two biological replicates and four technical-replicate libraries. Consequently, computational resampling of this matrix does not constitute independent biological replication.

### 5.2 Data Extraction and Reconstruction

Sparse contact records were extracted from the source `.hic` matrix and reconstructed into a midpoint-assigned contact representation.

The resulting dataset contains:

- 45,831 unique sparse contact records.
- 150,149 total observed contact counts.
- Contacts restricted to a maximum genomic separation of 32 kb.
- A corrected off-diagonal W80 profile comprising 122 genomic windows.

The analysis excludes the main diagonal and assigns contact pairs to genomic windows according to their midpoints.

Source records and experiment metadata:

- [ENCODE experiment ENCSR477GZK](https://www.encodeproject.org/experiments/ENCSR477GZK/)
- [ENCODE contact matrix ENCFF573OPJ](https://www.encodeproject.org/files/ENCFF573OPJ/)

### 5.3 Scope of Interpretation

The HCT116 analysis constitutes an initial real-data proof of concept.

It does not establish generalization across cell types, biological replicates, experimental protocols, or independent datasets.

All candidate regions identified from this dataset are treated as exploratory observations requiring further validation.

## 6. Exploratory Narrowing Candidates

### 6.1 Candidate Identification

Sliding-window analysis of the corrected HCT116 contact profile identified two candidate regions exhibiting locally reduced W80 values relative to comparison windows.

The initial detector used a descriptive narrowing-score threshold of 0.20. This threshold has not been calibrated to control statistical false discoveries.

### 6.2 Detected Genomic Regions

| Candidate | Genomic interval (GRCh38) | Peak W80 | Narrowing score |
| --- | --- | ---: | ---: |
| Region A | chr10:20,296,000–20,336,000 | 7.0 kb | 0.280 |
| Region B | chr10:20,720,000–20,792,000 | 6.5 kb | 0.387 |

The strongest individual analysis windows were:

- **Region A:** chr10:20,296,000–20,328,000.
- **Region B:** chr10:20,736,000–20,768,000.

Both peak scores were calculated against an initial reference W80 of 8.5 kb.

### 6.3 Scientific Interpretation

These results establish that the computational detector can identify localized deviations in observed Hi-C contact-distance profiles.

However, detection alone does not demonstrate changes in physical chromatin compaction, transcriptional regulation, or three-dimensional genome organization.

Subsequent analyses identified technical sensitivity associated with genomic mappability and contact visibility, particularly within Region B.

Therefore, both candidates remain exploratory and require additional validation before biological interpretation.

## 7. Statistical Robustness and Null-Model Evaluation

### 7.1 Contact-Thinning Experiments

To assess sensitivity to sequencing depth, observed contact counts were independently thinned using binomial sampling.

Each tested retention level was evaluated across 100 simulated realizations.

| Retained contact depth | Region A persistence | Region B persistence |
| --- | ---: | ---: |
| 25% | 66% | 97% |
| 10% | 55% | 94% |

Region B exhibited greater candidate-call persistence under contact thinning.

However, additional narrowing calls emerged at reduced depth. Candidate persistence therefore does not independently establish detector specificity or biological reproducibility.

### 7.2 Exploratory Null Models

The observed maximum narrowing score was compared against simulated maxima from several null-model configurations.

| Null-model configuration | Exceedances in 100 simulations |
| --- | ---: |
| Homogeneous distance-distribution null | 0 |
| Block-stratified null (100 kb) | 5 |
| Spatially smoothed null (50 kb) | 1 |
| Spatially smoothed null (100 kb) | 0 |
| Spatially smoothed null (200 kb) | 0 |

The null models produced different exceedance frequencies, demonstrating sensitivity to assumptions about spatial variation in contact intensity.

These simulations are exploratory and model-dependent. Their exceedance counts are not presented as calibrated p-values, genome-wide significance estimates, or formal false-discovery control.

### 7.3 Implications

The experiments demonstrate that apparent robustness under contact thinning is insufficient to establish that a detected narrowing event represents biological structure.

Reliable inference requires realistic null models, explicit data-quality criteria, appropriate multiple-testing procedures, and independent validation.

## 8. Genomic Mappability and Technical Artifacts

### 8.1 Independent Umap Analysis

To investigate whether apparent diagonal narrowing was influenced by sequence mappability, independent GRCh38 Umap tracks were examined at 50 bp and 100 bp read-length assumptions.

Region B exhibited reduced sequence mappability around several genomic bins with exceptionally low observed Hi-C contact support.

We subsequently recalculated W80 after excluding bins with low independently estimated mappability.

### 8.2 Mappability-Aware Results

| Analysis condition | Region B W80 | Reference W80 | Narrowing score |
| --- | ---: | ---: | ---: |
| Unmasked | 6.5 kb | 8.5 kb | 0.387 |
| Umap k50, threshold 0.50 | 7.5 kb | 8.5 kb | 0.181 |

The mappability-aware analysis reduced Region B's descriptive narrowing score substantially.

Region A also showed sensitivity to mappability masking, with its score decreasing from approximately 0.280 to 0.093 under the same k50 threshold.

### 8.3 Scientific Implications

These results establish that the apparent strength of candidate narrowing depends on the treatment of poorly mappable genomic positions.

However, the exact sequencing read length of the original experiment remains unverified. The Umap tracks are therefore treated as sensitivity analyses rather than definitive processing-matched filters.

The observed changes do not prove that all narrowing is artifactual. They demonstrate that sequence mappability is an important potential confounder requiring explicit treatment.

Reference: Karimzadeh et al. (2018), Umap and Bismap: quantifying genome and methylome mappability, Nucleic Acids Research.
https://doi.org/10.1093/nar/gky677

## 9. Mapping-Quality Filtering Sensitivity

### 9.1 Experimental Comparison

Two contact matrices from the same HCT116 experiment were compared:

- `ENCFF573OPJ`: mapping-quality-thresholded matrix (`inter_30`).
- `ENCFF726YRX`: more permissively filtered matrix (`inter_1`).

Both matrices represent pooled biological replicates from ENCODE experiment ENCSR477GZK.

### 9.2 Results

For the Region B peak window (chr10:20,736,000–20,768,000):

| Measurement | Stricter filtering | Relaxed filtering |
| --- | ---: | ---: |
| Total contact signal | 4,695 | 4,896 |
| Observed sparse contact pairs | 1,244 | 1,329 |
| Raw count-distribution W80 | 6.5 kb | 6.5 kb |
| Zero-support bins among eight examined | 6 | 0 |

Relaxed filtering increased contact visibility at all eight investigated low-support bins, while the targeted raw count-distribution W80 remained unchanged.

### 9.3 Interpretation and Limitations

These results indicate that local contact visibility is sensitive to mapping-quality filtering, even when the measured width remains stable.

However, this targeted comparison used raw contact counts summed by genomic separation. The primary W80 methodology uses mean contact intensity per eligible pair. These estimators are not generally equivalent.

An estimator-consistent comparison across candidate and matched-reference windows is therefore required before establishing robustness of the primary narrowing score.

This experiment measures technical processing sensitivity, not independent biological reproducibility.

## 10. Scientific Limitations

Although the initial computational experiments demonstrate promising methodological results, several limitations remain.

1. **Limited genomic scope:** Real-data validation has concentrated on a 1 Mb region of HCT116 chromosome 10. Generalization across the genome, cell types, and experimental protocols remains untested.

2. **Biological reproducibility:** The analyzed contact matrices combine biological replicates. Independent replicate-level validation has not yet been performed.

3. **Statistical calibration:** Candidate thresholds and exploratory null models do not yet provide calibrated genome-wide significance or false-discovery control.

4. **Sequence mappability:** The effective sequencing read length remains unverified, preventing definitive selection of a processing-matched mappability mask.

5. **Measurement reliability:** Sparse contact observations, low-visibility bins, and mapping-quality filtering can influence estimated widths and narrowing scores.

6. **Estimator consistency:** All sensitivity analyses must use equivalent W80 definitions and valid-pair normalization before their results can be compared directly.

7. **Biological interpretation:** No evidence currently establishes that the detected candidates represent physical chromatin compaction, regulatory boundaries, or causal biological mechanisms.

8. **Independent validation:** Additional genomic regions, biological datasets, and realistic technical controls are required.

These limitations motivate the next phase of method development and must be considered when interpreting the reported findings.

## 11. Reproducibility and Data Availability

### 11.1 Source Code

The repository contains the Python implementation of the local width estimator, candidate detector, quality-control procedures, synthetic benchmarks, and experimental analysis scripts.

Python dependencies are specified in `pyproject.toml` and `uv.lock`.

### 11.2 Experimental Data

The initial real-data experiments use publicly available HCT116 intact Hi-C data from ENCODE experiment ENCSR477GZK.

Original `.hic` files are not included in this repository because of their size. Their ENCODE accessions are provided in Section 5.

### 11.3 Computational Environment

Computationally intensive experiments were executed using Modal, with intermediate matrices and results stored in a persistent Modal Volume.

The repository contains scripts for these experiments, but a complete one-command reproduction workflow has not yet been established.

### 11.4 Research Artifacts

Generated artifacts include:

- Sparse HCT116 contact matrices.
- Corrected local W80 profiles.
- Candidate detection results.
- Coverage and contact-support diagnostics.
- Synthetic and null-model validation outputs.
- Genomic mappability sensitivity results.

These artifacts require a versioned public export and integrity manifest before the analysis can be considered independently reproducible.

### 11.5 Reproducibility Status

The core research code is version-controlled locally. Public GitHub synchronization, automated testing, and archival of the complete experimental results remain outstanding.

Accordingly, the current repository represents a documented research prototype rather than a fully reproducible scientific release.

## 12. Research Roadmap

The project is organized into four proposed research stages. These stages describe the broader research program rather than completed milestones.

### Stage 1 — Method Development and Proof of Concept

**Status: Substantial progress; not yet complete.**

Objectives:

- Implement and numerically validate the W80 estimator.
- Develop local-width scanning and candidate detection.
- Establish synthetic controls and contact-support diagnostics.
- Investigate exploratory narrowing candidates in real Hi-C data.
- Evaluate sensitivity to coverage, mappability, and mapping-quality filtering.
- Complete reproducible figures, automated tests, and research documentation.

### Stage 2 — Reliability and Statistical Calibration

**Status: Planned.**

Develop reliability-aware inference using explicit mappability masks, improved reference selection, realistic null models, calibrated detection thresholds, and uncertainty estimates.

### Stage 3 — Independent Biological Validation

**Status: Planned.**

Evaluate the framework across independent Hi-C datasets, biological replicates, genomic regions, cell types, and experimental conditions.

Investigate relationships between reliable contact-width changes and independently measured chromatin features without assuming causality.

### Stage 4 — Scientific Consolidation and Release

**Status: Planned.**

Integrate the validated methodology into a reproducible scientific software package.

Prepare comprehensive benchmarks, versioned datasets, publication-quality figures, technical documentation, and a research manuscript.

### Long-Term Research Objective

The broader objective is to develop a statistically defensible framework for identifying localized changes in Hi-C contact-distance structure while quantifying the effects of measurement uncertainty and technical artifacts.

Progression between stages depends on successful validation rather than a predetermined biological conclusion.

## 13. Authorship and Citation

**Researcher and developer:** Grold Otieno

Hi-C Diagonal Narrowing is an ongoing computational genomics research project.

The software and experimental results are currently under methodological development. No peer-reviewed publication, registered DOI, or formal software release is claimed at this stage.

Researchers using the underlying Hi-C data should acknowledge the ENCODE Consortium and cite the original data-generating studies where applicable.

Relevant resources:

- ENCODE experiment: https://www.encodeproject.org/experiments/ENCSR477GZK/
- Umap and Bismap: Karimzadeh et al. (2018), Nucleic Acids Research. https://doi.org/10.1093/nar/gky677

## 14. Scientific Disclaimer

The W80 statistic characterizes the genomic-distance distribution of observed Hi-C contacts. It does not directly measure three-dimensional chromatin compaction.

Detected narrowing candidates must not be interpreted as functional regulatory boundaries, structural transitions, or causal biological events without independent supporting evidence.

All reported candidate regions and statistical comparisons should be interpreted within the stated experimental assumptions and limitations.

**No biological significance is assumed from detection alone.**
