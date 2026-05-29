# Chapter: Introduction to Bayesian reasoning

## Overview

This chapter introduces the Bayesian modelling framework by contrasting it with the frequentist/maximum likelihood perspective. It develops the key ideas — Bayes' theorem, prior distributions, the posterior update, and Monte Carlo simulation — through a sequence of worked examples: significance testing, Covid-19 diagnostic testing, prior elicitation for a drug effectiveness parameter, convergence of expert opinion, and prior predictive (forward) sampling.

## Scripts

| File | Description |
|---|---|
| `bayesian-reasoning.R` | All R code for the chapter, in section order. |

## Code structure

### Significance testing (frequentist perspective)

Generates two simulated groups from Normal distributions (`n0=80`, `n1=78`) and carries out a two-sample t-test manually, computing the test statistic, degrees of freedom and one-sided p-value via `pt()`. The sampling distribution of T under H0 is plotted as a density curve with the p-value region shaded. This serves as the deductive-inference foil against which the Bayesian inductive approach is introduced.

### Bayes' theorem — Covid-19 testing

Applies Bayes' theorem to a diagnostic testing problem: given a test with 96% sensitivity and 95% specificity, the posterior probability of disease given a negative result is plotted as a function of the prior prevalence θ over the range [0, 1]. The key point is that the p-value analogy (falsely equating the sensitivity with the probability of disease given a negative result) ignores the prior and leads to the wrong conclusion when prevalence is low.

### Two experts

Demonstrates convergence of posterior opinions using a conjugate Beta-Binomial setup. Two experts start with Beta(3, 12) (low effectiveness) and Beta(12, 3) (high effectiveness) priors. After a small dataset (n=10) the posteriors remain distinct; after a large dataset (n=200) they are essentially identical. The visualisation uses a horizontal interval plot with arrows showing how each expert's 95% interval shifts at each update.

### Monte Carlo simulation

Shows how MC quantile estimates for a Normal(0, 1) distribution converge to the analytic values as sample size S grows from 10 to 1,000,000. The convergence is plotted on a log-scale x-axis, with horizontal dashed lines at the exact values `qnorm(0.025)` and `qnorm(0.975)`.

### Natural- vs original-scale parameters (Gamma example)

Illustrates prior elicitation on the natural (interpretable) scale when the model uses original-scale parameters. For Y ~ Gamma(nu, gamma), the natural-scale parameters are the mean `mu` and standard deviation `sigma`. Priors are placed on the natural scale (`mu ~ log-Normal(5.2, 0.2)`, `sigma ~ Exponential(0.35)`) and the relationships `gamma = sqrt(mu/sigma^2)` and `nu = mu*gamma` are applied to 10,000 MC draws to derive the implied priors for the original-scale parameters. Four density/density plots are produced for mu, sigma, nu and gamma.

### Drug example — Beta prior

The running "drug" example models the probability of back pain relief, theta, using a Beta(9.2, 13.8) prior. The parameters are derived by moment-matching: the prior information "most mass between 0.2 and 0.6" implies `mu=0.4`, `sigma=0.1`, which via the Beta moment equations gives `a=9.2`, `b=13.8`. The prior is simulated using `rbeta()` and summarised with `bmhe::stats()`. A tail-area probability `Pr(theta > 0.5)` is computed from the MC sample. The `bmhe` package installation commands are included.

### Prior information vs prior distribution (Beta vs logit-Normal)

Shows that the same prior information can be encoded in two equivalent ways: Beta(9.2, 13.8) directly on theta, or Normal(-0.405, 0.4137) on the logit scale. The logit-Normal density is implemented as a custom function `dlogitnorm()` using the change-of-variable formula. The near-equality of the two densities is verified numerically at theta=0.4 and visually by overlaying both density curves. `bmhe::logit()` and `bmhe::ilogit()` are used for the logit and inverse-logit transformations.

### Forward sampling (prior predictive)

Propagates the prior uncertainty on theta through a Binomial sampling model to generate the prior predictive distribution for y (number of successes in a hypothetical n=20 trial), using `rbinom()`. A threshold indicator `P.crit = (y >= 15)` is defined. The prior predictive is visualised as a bar chart with threshold-exceeding bars highlighted, and `bmhe::stats()` provides a numerical summary of `(P.crit, y, theta)` jointly.

## Dependencies

### R packages

```r
library(tidyverse)   # ggplot2 and data manipulation
library(bmhe)        # stats(), logit(), ilogit()
```

`bmhe` is the companion package for this book:

```r
# install.packages("remotes")
remotes::install_github("giabaio/bmhe_utils")
# Recommended alternative:
install.packages(
  "bmhe",
  repos = c("https://giabaio.r-universe.dev", "https://cloud.r-project.org")
)
```

### Data

All data are generated inline via `rnorm()`, `rbinom()`, `rbeta()`, etc. No external files are required.

### LaTeX / tikz note

The DAG figures (DGP, deductive/inductive inference, Covid testing decision tree) and several density plots use `dev: "tikz"` for typeset maths in axis labels. The extracted script uses the default ggplot2 device; axis labels render as plain strings rather than typeset maths.
