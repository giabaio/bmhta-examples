# Chapter: Survival analysis in HTA

## Overview

This chapter develops Bayesian survival modelling specifically for the HTA context, where extrapolation beyond observed trial follow-up is the primary objective. The running example is the digitised NICE TA174 dataset (RFC vs FC for chronic lymphocytic leukaemia), obtained from the `survHE` package. The chapter progresses from standard Kaplan-Meier estimation, through parametric Bayesian models fitted via `survHE`/Stan, to flexible M-spline models via `survextrap`, including incorporation of external aggregated data for extrapolation anchoring.

## Scripts

| File | Description |
|---|---|
| `survival.R` | All R code for the chapter, in section order. |

## Data

The TA174 dataset is bundled with `survHE` and loaded via `data(TA174, package="survHE")`. No external files are required for the main analyses. The chapter also uses a separate fictional dataset (`dat`) loaded from `data/survival/data_surv.Rdata` for some illustrative figures — this file is not publicly distributed.

The external aggregated data used in the M-spline external-data example are constructed inline as a `tibble` and represent hypothetical registry-style survival counts for the FC arm beyond the trial follow-up.

## Code structure

### Data preparation

The `TA174` dataset is restructured to a progression endpoint: individuals who die without progression, or are censored, are treated as censored (`status=0`); those who progress are events (`status=1`). The survival object is `Surv(time, status)` with `time` in months.

### Censoring illustration

A swim-lane diagram of five fictitious individuals illustrates right-censoring: filled circles = observed events, filled squares = study entry, and a vertical line = study close. Columns alongside show the observed time `t_i`, event indicator `d_i` and whether the true survival time is known.

### Kaplan-Meier curves

`survHE::fit.models()` with a Weibull distribution is used primarily to access the built-in KM computation and plotting facilities. The `plot()` method with `add.km=TRUE` overlays the KM step curves. Raw KM numbers (events, at risk, censored) are extracted from `m$misc$km` for tabulation.

### Generalised F model — Bayesian advantage

The Generalised F distribution (3 ancillary parameters: σ, p, q) illustrates why Bayesian inference is valuable under limited data. Under MLE, `p` is non-identifiable (95% CI spans essentially 0 to ∞). A Bayesian run with regularising priors (`log(p) ~ Normal(0, 0.5)`, `σ ~ Gamma(0.1,0.1)`, `q ~ Normal(0,2.5)`) produces a finite and clinically sensible estimate for `p`.

### survHE installation

The core `survHE` package installs from CRAN. The Bayesian HMC add-on (`survHEhmc`) and INLA add-on (`survHEinla`) install from GitHub or the r-universe repository. Only `survHEhmc` is required for the HMC analyses in this chapter.

### Fitting parametric Bayesian models

`fit.models()` with `method="hmc"` and `distr=c("wei","gom","lno")` fits Weibull (AFT), Gompertz and log-Normal models simultaneously using Stan via `rstan`. The resulting object `m` stores all three models and their DIC values. Key methods:
- `print(m, mod=k)` — posterior summaries in `survHE` notation (natural scale)
- `print(m, mod=k, original=TRUE, print_priors=TRUE)` — Stan parameter names, prior specs and convergence diagnostics (Rhat, n_eff)
- `rstan::traceplot(m$models[[k]])` — traceplot (direct rstan access)
- `rstan::stan_ac(m$models[[k]])` — autocorrelation plot

### Survival, hazard and cumulative hazard plots

`plot(m, ...)` returns a `ggplot` object. Key optional arguments:
- `add.km=TRUE` — overlay KM step curves
- `what="hazard"` or `what="cumhazard"` — hazard or cumulative hazard instead of survival
- `lab.profile=c("RFC","FC")` — custom treatment labels
- `mods=k` — restrict to one model
- `nsim=1000` — add a 95% uncertainty ribbon from the posterior (PSA)
- `t=seq(0,180)` — extrapolate to 180 months

When `nsim > 1` is used, `survHE` calls `make.surv()` internally to propagate the joint posterior uncertainty in all parameters to a distribution of survival curves.

Layer manipulation is demonstrated: `p$layers[[4]] = NULL` removes the CI ribbon around the KM curves, leaving only the KM step function alongside the parametric uncertainty band.

### Model selection

`model.fit.plot(m, type="DIC")` displays a bar chart of DIC values across models. Adding `scale="relative"` shows percentage increase from the best-fitting model, useful for communicating how much worse alternatives are.

### Separate modelling by treatment arm

When no single model fits both arms well, `fit.models()` is run separately on `data |> filter(treatment=="FC")` and `data |> filter(treatment=="RFC")`, with `~1` as the formula (intercept only). Both the location and ancillary parameters then differ between arms, effectively relaxing the PH assumption. `model.fit.plot(FC=m.ctl, RFC=m.trt, type="DIC", stacked=TRUE)` produces grouped bars for cross-arm comparison.

The combined extrapolation with `mods=c(1,5)` selects the 1st model from `m.ctl` (Weibull for FC) and the 5th overall = 2nd model from `m.trt` (Gompertz for RFC).

### M-spline models via survextrap

The `survextrap` workflow has three steps:

**1. `mspline_spec()`** — defines the spline structure: `df=6` basis functions, `add_knots=180` to extend the time grid to 180 months. Knot locations are placed at quantiles of observed event times.

**2. `survextrap()`** — fits the Bayesian model via Stan. The formula includes `treatment` as a covariate (PH assumption by default). Output includes: `alpha` (log baseline hazard scale φ), `coefs` (basis weights ω_k), `loghr` (log HR for treatment), `hr` (HR), `hsd` (smoothing SD σ). The smoothness prior σ ~ Gamma(2,1) penalises overfitting; as σ → 0, the hazard approaches constant.

**3. `survival()`** — extracts a tidy tibble of posterior survival summaries (`t`, `treatment`, `mean`, `2.5%`, `97.5%`). The M-spline curves are then combined with `survHE` plots using `geom_line()`, recoding the treatment column to match `survHE`'s internal `strata` naming convention.

### External aggregated data

`extdat` provides hypothetical Binomial counts `r_j` out of `n_j` individuals alive at `t_j^start`, for three future time intervals of the FC arm. Passing `external=extdat` to `survextrap()` anchors the extrapolation: the external data contribute to the likelihood through `pi_j = S(t_j^stop | theta) / S(t_j^start | theta)`. The resulting curves are pulled towards 0 at later times and show reduced uncertainty in the extrapolation region compared to the data-only M-spline.

## Dependencies

### R packages

```r
library(tidyverse)    # data manipulation and ggplot2
library(survHE)       # fit.models(), plot(), model.fit.plot(), make.surv(),
                      #   data(TA174), theme_survHE()
library(survextrap)   # mspline_spec(), survextrap(), survival()
library(rstan)        # traceplot(), stan_ac() (via survHEhmc in background)
```

Install the Bayesian survHE modules:

```r
install.packages("survHE")
install.packages(
  c("survHEhmc", "survHEinla"),
  repos = c("https://giabaio.r-universe.dev", "https://cloud.r-project.org")
)
install.packages("survextrap")
```

`survextrap` is available from CRAN and documented at [https://chjackson.github.io/survextrap/](https://chjackson.github.io/survextrap/).

### External software

Both `survHE` (HMC) and `survextrap` use Stan under the hood via `rstan`. Stan is installed automatically as a dependency but requires a C++ toolchain. See [https://mc-stan.org/rstan/](https://mc-stan.org/rstan/) for platform-specific setup instructions.

### Stan model code

The Stan code for the Weibull PH model (shown in the chapter) is embedded inside `survHEhmc`. Setting `save.stan=TRUE` in `fit.models()` saves the pre-compiled Stan code to the working directory, which can be used as a starting point for customisation.

### LaTeX / tikz note

Several figures use `dev: "tikz"` for axis labels with LaTeX maths. The extracted script uses the default ggplot2 device; axis labels render as plain strings rather than typeset maths.
