# Chapter: Aggregated level data and evidence synthesis

## Overview

This chapter develops Bayesian evidence synthesis using aggregated (summary-level) data from systematic literature review. The running example is the magnesium meta-analysis (S=16 RCTs), for which three models are compared: no-pooling, complete pooling and partial pooling (hierarchical random effects). A second worked example — neuraminidase inhibitors (NIs) for influenza prophylaxis — combines two independent evidence streams into a multiparameter evidence synthesis model and then carries the posterior uncertainty all the way to a cost-effectiveness analysis using BCEA.

## Scripts

| File | Description |
|---|---|
| `ald.R` | All R code for the chapter, in section order. |

## Code structure

### Setup

Global JAGS options suppress verbose output (`r2j.pb`, `r2j.quiet`, `r2j.print.program`). All three magnesium models are stored as named elements of a list `model` for easy comparison.

### Prior forward sampling

Before fitting any model, prior predictive checks verify that the chosen priors cover the full [0, 1] probability range on the implied scale:
- `alpha_s ~ Normal(0, sd=4)` implies a near-uniform marginal prior on `pi_s = logit^{-1}(alpha_s)`.
- Adding `delta_s ~ Normal(0, sd=2)` maintains the same coverage for the magnesium-arm probability.
- A demonstration of the Gamma(0.001, 0.001) precision prior shows why it is problematic for hierarchical SDs: almost all mass concentrates near tau=0, driving the SD towards implausibly large values.

### Magnesium meta-analysis — three models

All three models share the same Binomial sampling distribution (`r_st ~ Binomial(pi_st, n_st)`) and logistic linear predictor (`logit(pi_st) = alpha_s + delta_s*(Trt-1)`). They differ in how `delta_s` is modelled:

**No-pooling** (`model$no_pooling`). Each `delta_s` has an independent Normal(0, sd=2) prior. No pooled OR is estimated. Thinning (n.thin=4) compensates for slow mixing in small studies. Diagnostics flag nodes with n.eff < 400 (10% of nominal) and traceplots are shown for those.

**Complete pooling** (`model$complete_pooling`). All studies share a single log OR `d ~ Normal(0, sd=2)`. `OR = exp(d)` is monitored directly. This is the "fixed effects" model; it is dominated by ISIS-4 (the largest study) and cannot accommodate heterogeneity.

**Partial pooling** (`model$partial_pooling`). Study-specific log ORs `delta_s ~ Normal(d, sigma_delta)` are drawn from a common distribution. A PC prior `sigma_delta ~ Exponential(2.31)` encodes Pr(sigma_delta > 1) ≈ 0.1. Two pooled summaries are monitored: `OR = exp(d)` (pooled estimate) and `delta.pred ~ Normal(d, prec)` (predictive distribution for a new, exchangeable study). Convergence is substantially improved over no-pooling because information is shared across studies.

### DIC-based model comparison

`rjags::dic.samples()` is used to compute pD (the BUGS-style effective number of parameters) for all three models. The partial pooling model's pD is substantially lower than its nominal parameter count (35), quantifying the "borrowing of information" induced by the hierarchical structure. Results are presented in a summary table and a stacked bar chart decomposing DIC into mean deviance and penalty, for both pV and pD.

DIC-based model averaging weights are computed using `w_h = exp(-0.5 * ΔDIC_h) / Σ exp(-0.5 * ΔDIC_h)`. The partial pooling model dominates; the complete pooling model receives a modest non-zero weight because the DIC difference is in the ambiguous 5–7 range.

### Visualisations

**Forest plot.** Overlays posterior means and 95% intervals for `delta_s` from the no-pooling and partial pooling models alongside the complete pooling `d`, making the shrinkage visible for small studies.

**Shrinkage plot.** Scatter of `alpha_s` vs `delta_s` with arrows pointing from no-pooling to partial pooling estimates, labelled by trial number. Small studies show larger arrows (more shrinkage towards the grand mean `d`).

### Influenza — multiparameter evidence synthesis

Three versions of the model are fitted:

**Base-case** (`flu`). Alpha_s (head-to-head baselines) have independent Normal(0, sd=10) priors (no-pooling); `delta_s` (log ORs) are partially pooled via Normal(mu.delta, tau.delta); `beta_h` (single-arm incidence baselines) are partially pooled via Normal(mu.beta, tau.beta). The two modules are linked as `logit(p2) = logit(p1) + mu.delta`. PC priors `Exponential(1.65)` are used for all standard deviations (Pr(sigma > 1) ≈ 0.2).

**All-baselines-pooled** (`flu2`). Relaxes the assumption that head-to-head and single-arm baselines are independent by modelling all of them as `Normal(mu.beta, tau.beta)`. Appropriate if the two study sets are sufficiently similar.

**Baseline-risk adjusted** (`flu3`). Adds a term `gamma*(alpha_s - mean(alpha))` to the treatment arm linear predictor. `gamma ~ Normal(0, 0.01)` is monitored; a 95% interval spanning zero indicates no material baseline-risk effect in these data.

A coefficient plot from `bmhe::coefplot()` compares p1, p2 and OR across all three variants.

### Cost-effectiveness analysis

Posterior simulations from the base-case model are post-processed to compute population average effects and costs following the decision tree structure:

- Effects: `mu.e[t] = -l * p_t` (negative days of illness; NIs reduce illness duration).
- Costs: follow the two-branch decision tree (GP visit always; NIs cost added in treatment arm; influenza treatment cost added if infected).

`BCEA::bcea()` produces the economic summary at WTP = £1000/day. An empirical CDF overlay on the histogram of `Delta_e` allows direct reading of tail probabilities (e.g. Pr(Delta_e > 1 day)).

## Dependencies

### R packages

```r
library(tidyverse)    # data manipulation and ggplot2
library(R2jags)       # JAGS interface
library(BCEA)         # bcea(), summary.bcea()
library(bmhe)         # stats(), diagplot(), traceplot(), coefplot(), lognPar()
library(tinytable)    # tt(), style_tt(), theme_tt() for summary tables
```

`bmhe` is the companion package for this book:

```r
# install.packages("remotes")
remotes::install_github("<repo>/bmhe")
```

### External software

**JAGS** must be installed separately: [https://mcmc-jags.sourceforge.io](https://mcmc-jags.sourceforge.io)

### Data

All data are defined inline. There are no external data files required.

### LaTeX / tikz note

DAG figures (no-pooling, complete-pooling, partial-pooling, hierarchical, and the influenza decision tree) are rendered from `.tex` source files via the tikz engine. They have no corresponding R code and are not reproduced in this script.
