# Chapter: Indirect treatment comparisons (Network Meta-Analysis)

## Overview

This chapter extends Bayesian evidence synthesis to network meta-analysis (NMA), where indirect evidence from a connected network of studies can be combined with direct head-to-head evidence. The running example is the smoking cessation dataset (S=24 studies, T=4 interventions), used to illustrate fixed-effect and random-effect NMA, consistency checking between direct and indirect evidence, heterogeneity assessment, and a full cost-effectiveness analysis combining life-years gained with cost data.

## Scripts

| File | Description |
|---|---|
| `nma.R` | All R code for the chapter, in section order. |

## Data

The script requires `smoke.Rdata`, a list object containing:
- `r`: S × T_max matrix of quitter counts (NA where arm not present)
- `n`: S × T_max matrix of sample sizes
- `t`: S × T_max matrix of treatment indices
- `NT`: total number of treatment types

Update the `load()` path at the top of the script to point to the file on your system.

## Code structure

### Data wrangling

The wide-format data (one row per study) are pivoted to a long format where each row represents one study-arm. The variable `b` (baseline arm for each study) is computed as the minimum treatment index within each study. This long format is required for the nested-index JAGS model code (`delta[s[i], t[i]]`, etc.) and handles the unbalanced structure of the network (not all studies include all arms).

### Prior checks

A forward-sampling check verifies that `d[t] ~ Normal(0, sd=1)` for the log ORs implies a wide enough range on the OR scale (ORs easily exceeding 7), confirming the prior is not unduly restrictive.

### Helper functions

Three model functions and a data-filtering helper are defined before the main model runs:

**`direct_evidence()`**: a fixed-effect model restricted to the subset of studies directly comparing a specified pair (Active vs Baseline). Returns a pooled OR from direct evidence only, used for the consistency check.

**`heterogeneity()`**: a no-pooling model with study- and treatment-specific log ORs `d[s,t]`, producing one OR estimate per study per comparison. Used to visualise raw between-study heterogeneity.

**`make_data_no_pooling(type)`**: filters the long-format data to a specific pairwise comparison and runs either `direct_evidence()` or `heterogeneity()` across all six possible pairs. Returns a tidy tibble of posterior summaries.

### Fixed-effect NMA (`m_fe`)

Assumes a single pooled log OR `d[t]` per treatment arm (complete pooling across studies). Study-specific incremental effects `delta[s,t] = d[t] - d[b[s]]` are purely deterministic. All pairwise ORs are monitored in both directions. The absolute probability of quitting under each intervention is derived by linking `d[t]` to an informative prior on the baseline log-odds `rho ~ Normal(-2.6, precision=6.925)`, encoding an external estimate of 7% quitting without intervention (upper bound ~14%).

### Consistency check: direct vs indirect evidence

`make_data_no_pooling("direct")` is called to obtain direct-evidence-only OR estimates for all six pairwise comparisons. These are overlaid with the fixed-effect NMA estimates on a forest plot. Comparisons with large direct samples (C vs A: N=12846) are stable; pairs with sparse direct evidence (C vs B: N=255) show discrepancies, indicating inconsistency that motivates the random-effects model.

### Heterogeneity assessment

`make_data_no_pooling("heterogeneity")` produces study-specific ORs, visualised in a faceted forest plot with a log-scale x-axis. The large spread of estimates within each pairwise comparison confirms substantial between-study heterogeneity.

### Random-effects NMA (`m_re`)

Extends the fixed-effect model by making the study-specific incremental effects stochastic: `delta[s,t] ~ Normal(mu[s,t], tau)` with `mu[s,t] = d[t] - d[b[s]]`. A common standard deviation `sigma` is assumed across all treatment comparisons. The prior for `sigma` is a Half-Cauchy(0, 2), implemented in JAGS as a truncated t distribution:

```r
sigma ~ dt(0, 0.25, 1) %_% T(0, )
```

The `%_%` operator is a workaround required when passing JAGS model code as an R function — it is parsed to `T(0,)` by `R2WinBUGS::write.model()`. The script includes a helper function `ttrunc()` that computes the mean-precision t density in R (using the change-of-variable from R's standardised parameterisation) and a comparison plot of the Half-Cauchy vs an Exponential(0.4) PC prior, confirming they are closely aligned.

DIC comparison shows the random-effects model substantially outperforms the fixed-effect model, consistent with the observed heterogeneity. `rjags::dic.samples()` computes pD, the effective number of parameters, quantifying how much borrowing the hierarchical structure induces.

### Cost-effectiveness analysis

**Clinical benefit.** Life-years gained (LYG) from quitting smoking are modelled using external evidence from Mamun et al. (2004): `l_m ~ Normal(8.66, 0.495)` for men and `l_f ~ Normal(7.59, 0.679)` for women. Standard deviations are derived from the reported 95% intervals. Gender weights (0.567 men, 0.433 women) come from ONS 2024 UK smoking statistics. The intervention-specific expected benefit is `e[,t] = p[t] * l`, where `p[t]` comes directly from the random-effects NMA posterior.

**Costs.** Point estimates are based on NRT costs (35 patches × £1.30) plus clinic visit costs. Uncertainty is modelled as Gamma distributions with CV=0.3 (sd = 30% of mean), parameterised using `bmhe::gammaPar()`. No-intervention has zero cost.

**BCEA analysis.** `bcea()` is called with Group counselling (D) as reference. `multi.ce()` and `ceac.plot()` / `ceef.plot()` produce the probability-of-being-most-cost-effective curves and the efficiency frontier. The frontier shows Self-help is dominated for all values of the WTP threshold.

## Dependencies

### R packages

```r
library(tidyverse)      # data manipulation, ggplot2, pivot_longer
library(R2jags)         # JAGS interface
library(BCEA)           # bcea(), multi.ce(), ceac.plot(), ceef.plot()
library(bmhe)           # stats(), diagplot(), coefplot(), gammaPar(), lognPar()
library(directlabels)   # geom_dl() for inline curve labelling
library(scales)         # comma() for axis formatting
```

`bmhe` is the companion package for this book:

```r
# install.packages("remotes")
remotes::install_github("<repo>/bmhe")
```

### External software

**JAGS** must be installed separately: [https://mcmc-jags.sourceforge.io](https://mcmc-jags.sourceforge.io)

### LaTeX / tikz note

Several figures in this chapter use `dev: "tikz"` for axis labels with LaTeX maths. The extracted script uses default ggplot2 device; axis labels will render as plain strings rather than typeset maths.
