# Chapter: Markov models

## Overview

This chapter develops Bayesian cohort Markov modelling for HTA through two main examples. The first is a four-state HIV model (Chancellor et al. 1997) that combines Multinomial-Dirichlet transition probability estimation, partial-pooling evidence synthesis for the relative risk, a Weibull survival model for treatment-effect waning, matrix-algebra Markov simulation and discounted cost/LYG analysis. The second is a three-state cancer Markov model using the NICE TA174 digitised CLL data, linked to Bayesian Gompertz survival models via `survHE`.

## Scripts

| File | Description |
|---|---|
| `markov-models.R` | All R code for the chapter, in section order. |

## Data

**HIV example.** All data are defined inline as matrices and vectors. No external files are required.

**Three-state cancer model.** The TA174 dataset is loaded from `data/survival/ta174.rds`. This is a digitised and post-processed version of the NICE TA174 data (as described in the survival chapter). The path may need updating to match your local repository structure.

**PSM figure.** The partition survival analysis illustration loads pre-saved PFS/OS curves from `~/Dropbox/.../Pembro_PFS_OS.Rdata` (a hardcoded absolute path in the source). This figure has no standalone R computation beyond basic ggplot plotting and the relevant code is included for reference only.

## Code structure

### Dirichlet simplex visualisations

A custom `geom_simplex_canvas()` geom (adapted from `ggsimplex`) draws the equilateral triangle canvas with vertex labels a1, a2, a3. Six panels show different Dirichlet parameterisations: from flat/vague (a=0.85,0.85,0.85) through symmetric-concentrated (a=50,50,50) and various asymmetric configurations that pull mass towards specific vertices. Uses `brms::ddirichlet()` for the density computation and `ggsimplex::stat_simplex_density()` for the rendering.

### HIV example — data

The observed ZDV transition count matrix (4×4) is defined as `y`; inadmissible transitions (backwards progressions) are set to 0. Row-wise sample sizes `n` are derived with `apply(y, 1, sum)`. The four published RR studies are entered as point estimates and 95% CIs; log-scale values `x` and study SDs `sd` are derived analytically.

### Cost parameter elicitation

`bmhe::gammaPar(mean, sd)` is used to find Gamma(shape, rate) parameters matching the base-case and scenario cost estimates for direct care in each health state. Community care uncertainty is modelled by setting sd = 20% of the mean. Three `stat_function()` plots overlay the Gamma density, the 95% interval shaded region and the base-case to scenario range.

### JAGS model

The model function has four modules:
- **Multinomial-Dirichlet**: conjugate update of the ZDV transition matrix from count data.
- **Partial-pooling evidence synthesis**: random-effects Normal model on the log-RR across L=4 studies; `rho = exp(mu.theta)` is the pooled RR.
- **Weibull survival for treatment waning**: Binomial pseudo-counts `r[i]` out of `nr[i]` over time intervals encode expert assumptions; hazard parameterisation uses interval endpoints `u[i]`; `surv[j]` and `rho.star[j]` are derived for each Markov cycle.
- **Cost models**: Gamma priors for direct and community care per state; combination therapy drug cost `c.comb[j]` switches from `c_mono + c_3tc` to `c_mono` once `surv[j] < eps`.

### Post-processing — transition matrix

The time-varying 3TC+ZDV transition matrix is built using the logit-scale formulation (`eq-lambda2mat2`) for numerical stability: `logit(lambda^(2)) = logit(lambda^(1)) + log(RR) + adjustment`. A diagonal correction ensures row-sums equal 1. `which(lambda2 < 0, arr.ind=TRUE)` is used to detect any impossible transition probabilities.

### Markov model simulation

Arrays `m1` and `m2` (dimensions: n.sims × (J+1) × S) are filled via `m1[i,j,] = m1[i,j-1,] %*% lambda1[i,,]` and its equivalent for `m2` with the time-indexed `lambda2[i,,,j-1]`. The state occupancy trajectory is visualised as a faceted ribbon plot (mean ± 95% interval over MCMC draws).

### LYG and discounted costs

Occupancy probability arrays `pi1 = m1/1000` and `pi2 = m2/1000` are used to compute length-of-stay per state. LYG is the sum across non-Death states; the incremental effect `Delta_e = LYG2 - LYG1` is derived from the MCMC draws. Discounted costs use `mutate(across(starts_with("c_"), ~ .x/(1+d)^j))` to apply discounting in a single `dplyr` call. The cost tibbles are summarised using `group_by(sim) |> summarise(...)` to produce the incremental cost `Delta_c`.

### Three-state cancer Markov model

`survHE::make_data_multi_state(data)` restructures the TA174 event-history dataset into the long (mstate) format with one row per patient-transition. Three separate Gompertz models are fitted using `fit.models()` with `method="hmc"`, one for each competing risk set. `three_state_mm()` takes the three model objects and runs the Markov simulation for `nsim=1000` posterior draws over the time horizon `t=seq(0,120)`. `markov_trace()` produces the stacked area plot of state occupancy.

## Dependencies

### R packages

```r
library(tidyverse)    # data manipulation, ggplot2
library(R2jags)       # JAGS interface
library(survHE)       # fit.models(), make_data_multi_state(),
                      #   three_state_mm(), markov_trace()
library(bmhe)         # stats(), coefplot(), ilogit(), logit(), gammaPar()
library(ggsimplex)    # stat_simplex_density() for Dirichlet simplex plots
library(brms)         # ddirichlet() used inside ggsimplex
library(grid)         # for the custom simplex canvas geom
library(scales)       # scales::alpha() in the canvas geom
library(survHEhmc)    # fits models using survHEhmc
```

`bmhe` is the companion package:

```r
install.packages(
  "bmhe",
  repos = c("https://giabaio.r-universe.dev", "https://cloud.r-project.org")
)
```

`ggsimplex` installs from GitHub:

```r
remotes::install_github("marvinschmitt/ggsimplex")
```

### External software

**JAGS** must be installed separately: [https://mcmc-jags.sourceforge.io](https://mcmc-jags.sourceforge.io)

**Stan** is required by `survHEhmc` (called by `fit.models(..., method="hmc")`).

### LaTeX / tikz note

The Markov model structure DAGs (general Markov model, HIV DAG, three-state cancer DAG, partition survival DAG) are rendered from `.tex` source files via the tikz engine and have no corresponding R code.
