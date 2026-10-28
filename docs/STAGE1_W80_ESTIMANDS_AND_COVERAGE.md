# Stage 1: $W_{80}$ estimands and uncertainty validity

## Status

Methodological analysis supported by exact calculations,
automated tests, synthetic simulations, and exploratory
HCT116 contact data.

This document distinguishes computational validation,
conditional statistical guarantees, synthetic robustness,
and biological inference.

It does not claim validated biological confidence coverage.

## 1. Definition of $W_{80}$

Let $d=1,\ldots,D$ denote positive genomic separation
bands, with $K_d>0$ eligible genomic pairs at separation
$d$.

For observed aggregate contact count $Y_d$, define

```math
\widehat{\mu}_d=Y_d/K_d.
```


For a nonnegative intensity profile $\boldsymbol{\mu}$
with strictly positive total intensity, define

```math
W_{80}(\boldsymbol{\mu})
=
\min\left\{
w:
\frac{\sum_{d=1}^{w}\mu_d}
{\sum_{d=1}^{D}\mu_d}\geq 0.8
\right\}.
```


The empirical $W_{80}$ substitutes the distance-normalized
observed counts for the unknown intensities.

All definitions exclude diagonal contacts.

Midpoint-based window ownership is half-open. Pair
eligibility requires both genomic endpoints to pass the
specified mask, when a mask is used.

## 2. Conditional Poisson estimand

Assume fixed exposure counts $K_d$ and independent
distance-band counts

```math
Y_d\sim\mathrm{Poisson}(K_d\mu_d).
```


Here $\mu_d$ is the fixed underlying mean count per
eligible pair, and the target is

```math
W_{\mathrm{fixed}}=W_{80}(\boldsymbol{\mu}).
```


Conditional on total count $N=\sum_d Y_d>0$, the
distance-band counts follow a multinomial distribution with
probabilities proportional to $K_d\mu_d$.

The cutoff-specific exposure-bounded confidence procedure
in `poisson_uncertainty.py` uses this conditional structure,
with distance-specific exposure bounds and binomial
inversion.

Its confidence guarantee requires the model assumptions
and the specified fixed target. See
`docs/CUTOFF_EXPOSURE_W80_VALIDITY.md` for the formal
coverage argument and edge conditions.

Discreteness can make intervals conservative, meaning
actual coverage may exceed the nominal confidence level.

## 3. Fixed genomic-pair heterogeneity

Suppose pair $i$ at distance $d$ has a fixed intensity
$\lambda_{di}$, and pair counts are independent Poisson:

```math
X_{di}\sim\mathrm{Poisson}(\lambda_{di}).
```


Then

```math
Y_d=\sum_i X_{di}
\sim\mathrm{Poisson}
\left(\sum_i\lambda_{di}\right).
```


Consequently, unequal fixed pair-specific Poisson
intensities do not by themselves violate the
independent-Poisson band-sum model.

For this fixed map, the relevant intensity is

```math
\mu_d^{\mathrm{map}}
=
\frac{1}{K_d}\sum_i\lambda_{di}.
```


The associated target is the $W_{80}$ of this fixed map's
distance-normalized intensity profile.

Variation among genomic pairs is not automatically
biological overdispersion or evidence against Poisson
band-count sampling.

## 4. Random effects and distinct estimands

Consider an independently generated biological or
technical state $Z$ that changes the latent
distance-band intensities.

Conditional on $Z$, suppose

```math
Y_d\mid Z
\sim\mathrm{Poisson}(K_d\mu_d(Z))
```


independently across distance bands.

Two scientifically different targets arise.

**Realized-state target:**

```math
W_{\mathrm{latent}}(Z)
=
W_{80}\left(\boldsymbol{\mu}(Z)\right).
```


**Population-mean target:**

```math
W_{\mathrm{population}}
=
W_{80}\left(
\mathbb{E}_{Z}[\boldsymbol{\mu}(Z)]
\right).
```


In general,

```math
W_{80}\left(\mathbb{E}_{Z}[\boldsymbol{\mu}(Z)]\right)
\neq
W_{80}\left(\boldsymbol{\mu}(Z)\right).
```


The $W_{80}$ functional is nonlinear, and a random realized
state can move the target across discrete separation bins.

Even if conditional-Poisson intervals cover the latent
state-specific $W_{80}$, they need not cover the $W_{80}$ of the
population-mean profile.

These are different estimands, not interchangeable labels.

## 5. Shared versus distance-specific effects

A single positive multiplicative effect shared across
all distance bands preserves their relative intensities.

It therefore preserves the latent $W_{80}$, and conditioning
on the total number of contacts retains the relevant
multinomial probabilities.

By contrast, distance-specific or distance-block effects
change relative intensities and can shift latent $W_{80}$.

The unconditional distribution may also exhibit
dependence or overdispersion.

The consequences depend on both the stochastic mechanism
and the scientific target.

## 6. Synthetic stress-test evidence

### 6.1 Model misspecification

The experiment

`experiments/stress_test_w80_model_misspecification.py`

used 1,000 datasets per model and profile.

Selected empirical population-target coverage:

| Generating model | A-like | B-like |
| --- | ---: | ---: |
| Independent Poisson | 100.0% | 98.6% |
| Shared global Gamma depth | 99.8% | 98.8% |
| Independent pair Gamma-Poisson | 95.8% | 84.8% |
| Correlated eight-distance blocks | 41.0% | 34.2% |

These are Monte Carlo results for specified synthetic
mechanisms. They are not fitted HCT116 models.

### 6.2 Explicit target comparison

The experiment

`experiments/diagnose_w80_estimand_under_heterogeneity.py`

evaluated population and latent targets separately.

For Region B-like eight-distance-block random effects:

- Population-target coverage: 36.5%.
- Latent-target coverage: 99.3%.
- Fraction of simulations with shifted latent $W_{80}$: 79.8%.

This is consistent with target mismatch playing a major
role in the simulated population-coverage failure.

### 6.3 Independent fixed maps

The experiment

`experiments/audit_fixed_heterogeneity_maps.py`

generated 40 independent fixed intensity maps for each
synthetic profile, with 80 independent Poisson observations
from each map.

The unit of independent variation between intensity maps
is the map, not each Poisson observation.

A subsequent cluster bootstrap

`experiments/audit_w80_map_level_uncertainty.py`

resampled the 40 maps 10,000 times.

| Measure | A-like | B-like |
| --- | ---: | ---: |
| Fixed-map $W_{80}$ coverage | 99.6% | 99.5% |
| Reference population $W_{80}$ coverage | 97.2% | 83.6% |
| Coverage difference | 2.4 pp | 15.9 pp |
| 95% bootstrap interval for difference | 0.5–5.6 pp | 7.3–25.5 pp |
| Shifted maps | 8/40 | 15/40 |

The Region B-like population-target coverage had a
map-level bootstrap interval of 74.0–92.2%.

These bootstrap intervals quantify uncertainty under the
specified synthetic map-generating distribution; they do
not quantify uncertainty for the biological HCT116 system.

The map-level sample size remains only 40 per profile.

## 7. Existing real HCT116 evidence

The analyzed source is the pooled HCT116 GRCh38 Hi-C
file ENCFF573OPJ, chr10:20–21 Mb, 500-bp resolution.

The contact values are integer-valued despite float64
storage. Source checksums were verified.

All 122 corrected $W_{80}$ windows were reproduced from the
original sparse contact representation.

A descriptive real-data diagnostic reported median
within-distance pair Fano ratio approximately 1.09.
This is not a formal test of Poisson band-count sampling,
because pair-specific expected rates are unknown.

At full depth, the model-conditional $W_{80}$ intervals had
median span three 500-bp bins.

Depth-matched binomial thinning yielded median interval
spans of six bins at 25% retention and nine bins at 10%
retention. Inclusion of the original full-depth empirical
$W_{80}$ was approximately 99.2% and 98.5%, respectively.

Artificial thinning is conditional on an already observed
pooled matrix and cannot estimate between-biological-
replicate variation or establish biological 95% coverage.

The source experiment includes biological replicates 1
and 2, but the available metadata inventory did not
identify a released single-replicate processed Hi-C
matrix. Replicate-level BAM alignments exist.

Mappability-based bin exclusion has been explored, but
genomic missing-bin status is not definitively established.

## 8. Permissible and impermissible interpretations

Supported:

- $W_{80}$ and its genomic pair geometry are computationally
  reproducible under the tested implementations.
- The cutoff-specific interval has conditional guarantees
  under its stated statistical model.
- Synthetic random-effect mechanisms can cause serious
  failures when fixed-map intervals are interpreted as
  intervals for a population-mean $W_{80}$.
- Region B shows stronger technical candidate stability
  than Region A under the investigated thinning protocol.

Not supported:

- Claiming biological 95% confidence coverage for HCT116.
- Treating 122 overlapping windows as 122 biological
  replicates.
- Treating simulated intensity maps as observed
  biological replicates.
- Claiming significance for candidate regions based only
  on their individual $W_{80}$ intervals.
- Claiming pair-level Fano ratios identify a particular
  biological overdispersion model.
- Treating the $W_{80}$ of one pooled sample as automatically
  equal to the $W_{80}$ of a broader biological population.

## 9. Requirements for stronger biological inference

Further work should define the scientific target before
selecting or calibrating an interval.

For inference about an observed state's latent contact
profile, assumptions about contact sampling, pair
eligibility, dependence, and technical artifacts must be
examined.

For inference about population-level variation, independent
biological samples or a separately justified hierarchical
model are needed.

Candidate-versus-reference comparisons also require a
joint analysis of shared genomic contacts, dependence
across windows, and exploratory candidate selection.

Until these requirements are addressed, label intervals
as *model-conditional* and findings as exploratory.
