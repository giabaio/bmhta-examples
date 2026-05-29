# Chapter: Bayesian Software

## Overview

This chapter covers the practical workflow for running Bayesian models from R using JAGS and BUGS. It introduces the BUGS model syntax, explains how to interface R with JAGS via `R2jags`, walks through post-processing MCMC output (traceplots, autocorrelation, convergence diagnostics), documents the key default-setting differences between BUGS and JAGS that affect results even when nominal arguments are identical, and shows how to implement non-standard distributions in JAGS using the "zero-trick".

## Scripts

| File | Description |
|---|---|
| `bayesian-software.R` | All R code for the chapter, in section order. |

## External model files

JAGS model code can be passed either as a path to a plain-text `.txt` file or as an R function (the approach used throughout this script). The drug example model is shown as an R function (`model.code`) and the PC prior example as `model_pc_bern`. If you prefer external files, save the model blocks to `.txt` and update the `model.file` argument in the `jags()` call accordingly.

The eight-schools model used in the BUGS vs JAGS comparison is bundled with `R2OpenBUGS` and located automatically via `system.file()`.

## Code structure

**Package installation.** Installs `R2jags` and `R2OpenBUGS`. Only one is strictly needed; this book uses JAGS throughout.

**Drug example — JAGS run.** Fits a Beta-Binomial model for a drug efficacy pilot study. The prior is Beta(9.2, 13.8), the observed data are y=15 successes out of m=20 trials and the model also produces a posterior predictive distribution for y.pred (successes in a future trial of size n=40) and a binary indicator P.crit (whether y.pred meets a critical threshold of 25). Key `jags()` arguments are documented in-line.

**Post-processing.** Demonstrates the standard post-processing workflow using `bmhe` helper functions:
- `bmhe::traceplot()` for visual convergence inspection.
- `bmhe::acfplot()` for autocorrelation.
- `bmhe::diagplot()` for R-hat and effective sample size.
- `bmhe::posteriorplot()` for posterior density plots.

**BUGS vs JAGS defaults.** Runs the eight-schools hierarchical model with both `R2OpenBUGS::bugs()` and `R2jags::jags()` using identical top-level arguments, then compares the results. The key practical differences are:
- *Thinning*: BUGS defaults to `n.thin=1`; JAGS defaults to `n.thin=max(1, floor((n.iter-n.burnin)/1000))`. With `n.iter=5000` and `n.chains=3` this means BUGS retains 7500 simulations and JAGS retains 3750 by default.
- *Object structure*: `R2jags` wraps all output inside a nested `$BUGSoutput` slot; `R2OpenBUGS` places results at the top level.
- *Initial values*: BUGS draws from the prior; JAGS uses the prior mean/mode.
- *Random seed*: BUGS uses 1; JAGS uses 123.

**Zero-trick (PC prior).** Implements the penalised complexity (PC) prior for a Bernoulli parameter inside JAGS using the zero-trick: a pseudo-observation w=0 is modelled as Poisson(phi), where phi encodes the negative log of the target density. A large constant C=10000 keeps phi positive. The model is run with no observed data (y=NA) to verify that the implied marginal distribution for theta matches the analytically derived PC prior, and Pr(theta > 0.75) is estimated by Monte Carlo.

## Dependencies

### R packages

```r
library(tidyverse)    # ggplot2 and data manipulation
library(R2jags)       # R interface to JAGS
library(R2OpenBUGS)   # R interface to OpenBUGS (comparison only)
library(bmhe)         # book companion package: traceplot(), acfplot(),
                      #   diagplot(), posteriorplot(), stats()
```

`bmhe` is the companion package for this book:

```r
# install.packages("remotes")
remotes::install_github("<repo>/bmhe")
```

### External software

- **JAGS** must be installed separately: [https://mcmc-jags.sourceforge.io](https://mcmc-jags.sourceforge.io)
- **OpenBUGS** (optional, comparison section only): [https://www.mrc-bsu.cam.ac.uk/software](https://www.mrc-bsu.cam.ac.uk/software)

### LaTeX / tikz note

The DAG figures in this chapter (the BUGS DAG and the predictive distribution DAGs) are rendered from `.tex` source files via the `tikz` engine in the original document. These figures have no corresponding R code and are not reproduced in the script.
