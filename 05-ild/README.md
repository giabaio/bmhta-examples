# Chapter: Cost-effectiveness analysis with individual-level data

## Overview

This chapter develops a full Bayesian cost-effectiveness analysis pipeline for the 10TT trial (synthetic data), a two-arm RCT of a brief weight-loss intervention for obese adults. Starting from raw individual-level HRQL and resource-use data, the code builds three progressively more realistic joint models for the bivariate outcome (QALYs, total costs), compares them using DIC and WAIC, and feeds the posterior simulations into BCEA to produce cost-effectiveness summaries and a structural PSA via model averaging.

## Scripts

| File | Description |
|---|---|
| `ild.R` | All R code for the chapter, in section order. |

## Data

The script expects a CSV file at `data/ild/10TT_synth_280921.csv`. This is a synthetic dataset mimicking the original 10TT trial structure. Columns include individual ID, treatment arm, demographic variables (age, sex, BMI, GP visits), HRQL utility scores at baseline and 3, 6, 12, 18 and 24 months (`qol_0`, `qol_3`, ..., `qol_24`) and total costs.

## Code structure

### Data preparation

**QALY computation.** Pivots the wide-format utility data to long format, applies a 3.5% annual discount rate (NICE default) using a helper function `disc()`, and computes QALYs per individual as the discounted area under the piecewise-linear utility curve (trapezoid rule). Discounted QALYs are merged back into the main dataset.

**Data objects for JAGS.** Costs are scaled to £1000 units before passing to JAGS, to keep the range of data compatible with the Normal(0, 100) priors on regression coefficients. The baseline utility `u0star` is mean-centred. All objects are collected in a named list `data`.

### Statistical models

All three models are specified as R functions and passed to `R2jags::jags()`. Results are stored as named elements of a list `model` to keep everything in one place.

**Normal/Normal independent** (`model$nn_indep`). Separate Normal linear regressions for effects and costs, controlling for treatment arm and centred baseline utility. Arm-specific precision parameters. Serves as the baseline comparator. PC priors on standard deviations: Exponential(5.75) for effects (Pr(σ_e > 0.8) ≈ 0.01) and Exponential(0.35) for costs on the £1000 scale (Pr(σ_c > 2) ≈ 0.5).

**Normal/Normal MCF** (`model$nn_mcf`). Extends the independence model by regressing costs on (centred) effects within each arm: `phi.c[i] = beta0 + beta1*(Trt[i]-1) + beta2*(e[i]-mu.e[Trt[i]])`. Under joint Normality this is algebraically equivalent to the SUR model. The marginal cost SD and correlation coefficient can be derived analytically and are either monitored inside JAGS or reconstructed post-hoc in R using `lambda.c`, `sigma.e` and `beta2`.

**Gamma/Gamma MCF** (`model$gg_mcf`). Replaces Normal sampling distributions with Gamma, appropriate for the observed right-skewed costs and left-skewed QALYs. QALYs are first transformed as `estar = 3 - e` to flip the skew; the back-transformation `mu.e[t] = 3 - mustar.e[t]` is valid because the transformation is linear. Log-linear predictors ensure positive means. PC priors on the Gamma shape parameters: Exponential(0.15) (Pr(ν > 30) ≈ 0.01).

### Posterior predictive g-computation

Demonstrates the general posterior predictive approach for recovering population average outcomes when the transformation from the model scale to the natural scale is non-linear. For each MCMC draw, `N_mc = 4000` values of `estar` are simulated from the fitted Gamma distribution and back-transformed; the mean across these draws gives one posterior sample of the population average QALY. Equivalent to the direct JAGS calculation when the transformation is linear.

### Model selection

**Manual pD.** Reconstructs Dhat = D(E[θ]) by evaluating the Gamma log-likelihood at the posterior mean parameters for both the effects and costs components, then computes pD = Dbar - Dhat.

**pD via R2jags.** Re-runs all three models with `pD=TRUE`, which calls `rjags::dic.samples()` automatically after the main MCMC run.

**DIC comparison table.** Assembles a table of pV, DIC (pV-based), pD and DIC (pD-based) for the three models. The Gamma/Gamma MCF is far superior (ΔDIC > 100 relative to the Normal/Normal models).

**WAIC and LOO-CV.** A variant of `gg_mcf` adds individual log-likelihood contributions (`log.lik[i]`) using the JAGS built-in `logdensity.gamma()`. The resulting matrix is fed to `loo::waic()` and `loo::loo()`.

**DIC weights.** Implements the weighting formula w_h = exp(-0.5 * ΔDIC_h) / Σ exp(-0.5 * ΔDIC_h) and plots the decay of a model's weight as a function of its ΔDIC, showing that weights effectively vanish beyond ΔDIC ≈ 10.

### Cost-effectiveness analysis

Three `bcea()` objects are constructed from the MCMC simulations for population average effects and costs. The `graph="gg"` option returns ggplot objects, whose underlying data are extracted and stacked to produce a single overlaid CEAC plot for all three models. `contour2()` produces the cost-effectiveness plane with joint posterior contour.

`BCEA::struct.psa()` computes DIC-weighted model-averaged simulations for effects and costs and returns a `bcea` object, which can be used with all standard BCEA methods.

## Dependencies

### R packages

```r
library(tidyverse)    # data manipulation, ggplot2, pivot_longer
library(R2jags)       # JAGS interface
library(BCEA)         # bcea(), ceac.plot(), contour2(), struct.psa()
library(bmhe)         # stats(), coefplot(), lognPar(), gammaPar()
library(loo)          # waic(), loo()
library(tinytable)    # tt(), style_tt() for the DIC summary table
```

`bmhe` is the companion package for this book:

```r
# install.packages("remotes")
remotes::install_github("<repo>/bmhe")
```

### External software

**JAGS** must be installed separately: [https://mcmc-jags.sourceforge.io](https://mcmc-jags.sourceforge.io)

### LaTeX / tikz note

Several figures in this chapter use `dev: "tikz"` for axis labels with LaTeX maths. The extracted script drops the tikz device requirement; all code runs with the default ggplot2 device, but axis labels will appear as plain strings rather than typeset maths.
