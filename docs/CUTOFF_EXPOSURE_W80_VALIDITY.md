# Statistical validity of cutoff-specific $W_{80}$ confidence intervals

## 1. Population estimand

For contact separations $d=1,\ldots,D$, define:

- $K_d>0$: fixed eligible genomic-pair count.
- $\mu_d\geq0$: population mean contact intensity per eligible pair.
- $S=\sum_{d=1}^{D}\mu_d>0$.

The population cumulative intensity profile is

```math
F(d)=\frac{\sum_{j=1}^{d}\mu_j}{S}
```


The population $W_{80}$ is

```math
w=\min\{d:F(d)\geq 0.80\}
```


This is a contact-distance statistic, not a direct measurement
of three-dimensional chromatin compaction.

## 2. Sampling assumptions

Assume independent band totals:

```math
Y_d\sim\mathrm{Poisson}(K_d\mu_d)
```


Let $N=\sum_{d=1}^{D}Y_d$.

Conditional on $N=n>0$, the band counts are multinomial,
with probabilities proportional to $K_d\mu_d$.

Consequently, the observed cumulative count $C_d$ satisfies

```math
C_d\mid N=n\sim\mathrm{Binomial}(n,Q_d)
```


where

```math
Q_d=\frac{\sum_{j=1}^{d}K_j\mu_j}{\sum_{j=1}^{D}K_j\mu_j}
```


No equal-exposure assumption is necessary.

## 3. Cutoff-specific exposure bounds

For every cutoff $d<D$, partition the distances into:

- Prefix: separations $\leq d$.
- Tail: separations $>d$.

Let the minimum and maximum exposures within each
partition be $K_{\mathrm{prefix,min}}$, $K_{\mathrm{prefix,max}}$,
$K_{\mathrm{tail,min}}$, and $K_{\mathrm{tail,max}}$.

Set $t=0.80$.

If $F(d)\geq t$, then

```math
Q_d\geq p_{\mathrm{low}}(d)
```


where

```math
p_{\mathrm{low}}(d)=\frac{tK_{\mathrm{prefix,min}}}{tK_{\mathrm{prefix,min}}+(1-t)K_{\mathrm{tail,max}}}
```


If $F(d)\leq t$, then

```math
Q_d\leq p_{\mathrm{high}}(d)
```


where

```math
p_{\mathrm{high}}(d)=\frac{tK_{\mathrm{prefix,max}}}{tK_{\mathrm{prefix,max}}+(1-t)K_{\mathrm{tail,min}}}
```


These inequalities follow from bounding the
exposure-weighted prefix and tail intensities.

Neither probability must straddle 0.80.

## 4. Confidence interval construction

For confidence level $1-\alpha$ and observed $N=n>0$:

```math
a=\frac{\alpha}{2}
```


At each cutoff $d<D$, calculate

```math
\begin{aligned}
r_{\mathrm{low}}(d)&=\mathrm{BinomialPPF}(a;n,p_{\mathrm{low}}(d))\\
r_{\mathrm{high}}(d)&=\mathrm{BinomialISF}(a;n,p_{\mathrm{high}}(d))
\end{aligned}
```


Define $L$ as the earliest cutoff satisfying

```math
C_d\geq r_{\mathrm{low}}(d)
```


If none exists, set $L=D$.

Define $U$ as one plus the latest cutoff satisfying

```math
C_d\leq r_{\mathrm{high}}(d)
```


If none exists, set $U=1$.

The resulting interval is $[L,U]$.

The construction does not assume that the
cutoff-specific thresholds are monotonic.

## 5. Conditional coverage argument

Let $w$ be the true population $W_{80}$.

When $w<D$, $F(w)\geq t$ implies

```math
Q_w\geq p_{\mathrm{low}}(w)
```


Binomial stochastic ordering and the lower-tail
quantile definition imply

```math
\Pr(C_w<r_{\mathrm{low}}(w)\mid N=n)\leq\frac{\alpha}{2}
```


If $C_w\geq r_{\mathrm{low}}(w)$, the cutoff $w$ is an admissible
lower-bound candidate. Therefore $L\leq w$.

Thus

```math
\Pr(L>w\mid N=n)\leq\frac{\alpha}{2}
```


When $w>1$, $F(w-1)<t$ implies

```math
Q_{w-1}\leq p_{\mathrm{high}}(w-1)
```


The corresponding upper-tail bound implies

```math
\Pr(C_{w-1}>r_{\mathrm{high}}(w-1)\mid N=n)\leq\frac{\alpha}{2}
```


If $C_{w-1}\leq r_{\mathrm{high}}(w-1)$, the cutoff $w-1$ is
an admissible upper-bound candidate. Therefore $U\geq w$.

Thus

```math
\Pr(U<w\mid N=n)\leq\frac{\alpha}{2}
```


At the boundary values $w=1$ and $w=D$, the
respective one-sided error is impossible.

Applying the union bound yields

```math
\Pr(L\leq w\leq U\mid N=n)\geq1-\alpha
```


for every $n>0$ under the specified sampling model.

Only the two cutoffs adjoining the true $W_{80}$ are
needed for this argument, not simultaneous
testing across every separation.

## 6. Limitations

This guarantee depends on:

- Independent Poisson contact-band totals.
- Fixed, positive eligible-pair exposures.
- Raw nonnegative integer contact counts.
- A positive observed total contact count.
- A population $W_{80}$ defined by mean intensity.

The implementation abstains when no contacts are observed.

The guarantee does not automatically transfer to:

- Biological Hi-C with contact dependence or overdispersion.
- Data-dependent missing-bin masks.
- Normalized or transformed contact matrices.
- Selected genomic candidates or multiple-window inference.
- Comparisons of $W_{80}$ between two genomic windows.

Independent biological validation remains necessary.

## 7. Computational evidence

The reusable implementation matched the independent
experimental implementation on 7,200 synthetic trials,
with zero mismatches.

A separate exact-enumeration audit evaluated 144
small-model configurations. All achieved at least
95% conditional coverage.

These tests support correctness but do not replace
the statistical argument or validate its assumptions
for biological Hi-C.

## 8. Implementation and verification

Implementation:
`src/hic_diagonal_narrowing/poisson_uncertainty.py`

Function:
`cutoff_exposure_w80_interval_from_bands`

Tests:
`tests/test_cutoff_exposure.py`

Experiments:
`experiments/compare_cutoff_exposure_rank.py`
`experiments/verify_cutoff_exposure_rank_exact.py`
`experiments/verify_cutoff_exposure_implementation.py`
