# Chapter: Learning from data — Bayesian computation

## Overview

This chapter covers the computational machinery needed to go from a prior distribution and observed data to a posterior distribution. It develops the material in three stages: choosing prior distributions, computing the posterior analytically when structure allows it, and approximating it numerically via Markov Chain Monte Carlo (MCMC) when it does not. A running example throughout is a drug-effectiveness experiment modelled as Binomial with an unknown success probability θ.

## Scripts

| File | Description |
|---|---|
| `bayesian-computation.R` | All R code extracted from the chapter, in section order. |

## Code structure

The script follows the chapter's section order.

**Setup.** Defines the global number of Monte Carlo simulations (`nsim = 10000`) and the Beta(9.2, 13.8) prior for the running drug example. These objects are reused in several later sections.

**Vague / Uniform priors.** Compares the Beta(9.2, 13.8) and Uniform(0,1) priors through (i) density plots, (ii) forward-sampling (prior predictive) simulations, (iii) analytically computed Beta-Binomial prior predictive distributions, and (iv) the closed-form posterior after observing y=15 successes out of n=20. Monte Carlo summary uses `bmhe::stats()`.

**Prior elicitation via pseudo-data.** Constructs four Beta posteriors by "pretending" to have seen different (y₀, n₀) pseudo-datasets, starting from a Uniform(0,1). This gives a principled way to calibrate a prior to match a desired level of prior uncertainty.

**Vague Normal prior.** Plots a Normal(0, 100000) prior — proper but extremely diffuse — as an alternative vague prior for unbounded parameters (e.g. regression coefficients).

**Conjugate analysis — Beta-Binomial.** Plots the prior, (rescaled) likelihood and conjugate posterior for the drug example. The posterior is Beta(a₀+y, b₀+n−y); no simulation is needed. Uses `directlabels::geom_dl()` to annotate curves directly.

**Covid-19 vaccine example.** Replicates the Bayesian sample-size and efficacy analysis from the Pfizer/BioNTech Phase II/III trial. Vaccine efficacy is reparameterised as θ = π_vac/(π_vac + π_plac) to allow a single Beta prior. Observed arm counts are rescaled to restore the 1:1 allocation assumption, then the Beta prior is updated conjugately and 100,000 Monte Carlo draws characterise the posterior for VE.

**Penalised Complexity (PC) priors.** Implements the PC prior for two cases:
- *Binomial parameter*: `pc_prior()` computes the KLD-based PC prior for a Bernoulli θ relative to a base model at θ₀; numerical integration via `integrate()` calibrates the tail probability.
- *Standard deviation*: Derives an Exponential prior for σ and plots both the type-2 Gumbel prior on τ and the Exponential prior on σ under the constraint Pr(σ > 1) = 0.01.

**Jeffreys' prior (Binomial example).** Defines a custom legend-key function (`draw_key_numbered`) and a Binomial likelihood function (`lik()`), then overlays Jeffreys' Beta(0.5, 0.5) prior, the likelihood and the resulting posterior for four (y, n) combinations.

**Gibbs sampling — semi-conjugated Normal model.** Implements a Gibbs sampler from scratch for the model y ~ Normal(μ, σ²) with Normal prior on μ and Gamma prior on τ = 1/σ². The full conditionals are available analytically (Normal for μ, Gamma for τ). Key steps: data input, prior hyper-parameter specification, random initialisation, and the Gibbs loop. An animated step-by-step visualisation uses `mvtnorm::dmvnorm()` contours; `ggExtra::ggMarginal()` adds marginal histograms.

**MCMC diagnostics.** Loads pre-saved chain output (`data/bayesian-computation/convergence.RData`) and implements:
- Traceplot with vertical burn-in marker.
- `Rhat()`: Potential Scale Reduction (Brooks-Gelman-Rubin R̂). Values close to 1 indicate convergence.
- `n_eff()`: Effective Sample Size approximation. Values well below the nominal sample size signal autocorrelation.
- Monte Carlo Standard Error (MCSE = sd / √ESS).
- Autocorrelation function plots via `bmhe::acfplot()`.

The diagnostics are also applied to the semi-conjugated model run with two chains from different starting points.

**Hamiltonian Monte Carlo illustration.** Plots a schematic of HMC dynamics: the target density and its negative log-density ("potential energy") on a dual axis, and a discretised Hamiltonian trajectory from a starting position with high potential energy.

## Dependencies

### R packages

```r
library(tidyverse)     # core data manipulation and ggplot2
library(directlabels)  # direct curve labelling (geom_dl)
library(mvtnorm)       # multivariate Normal density for contour overlays
library(ggExtra)       # marginal histograms (ggMarginal)
library(bmhe)          # book companion package: stats(), acfplot()
```

`bmhe` is the companion package for this book. Install from GitHub if not already available:

```r
# install.packages("remotes")
remotes::install_github("<repo>/bmhe")
```

### External data

- `data/bayesian-computation/convergence.RData` — pre-run MCMC output (two chains) used for the convergence and autocorrelation diagnostic figures.

### LaTeX / tikz note

In the source `.qmd` file, most figures are rendered with `dev: "tikz"`, which passes the plot through a LaTeX engine to render maths in axis labels and annotations (e.g. `$\\theta$`, `$\\mathcal{L}$`). This requires a working LaTeX installation and the R package `tikzDevice`.

The extracted `.R` script drops the `tikz` device requirement: when run interactively or in a standard Quarto/knitr context, ggplot2 will render the labels as plain text strings rather than LaTeX. If you want matching typeset output, add `options(tikzLatex = "<path-to-latex>")` and use `tikzDevice::tikz()` as the graphics device.
