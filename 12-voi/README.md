# Chapter: Value of Information

## Overview

This chapter introduces Value of Information (VoI) analysis as a principled framework for HTA uncertainty analysis and research prioritisation. It covers the three main VoI measures — EVPI, EVPPI and EVSI — along with the Expected Net Benefit of Sampling (ENBS) for sample size determination. The chapter is mainly methodological; the code implements self-contained illustrations of each concept and uses built-in datasets from `BCEA` and `voi` for the applied examples.

## Scripts

| File | Description |
|---|---|
| `voi.R` | All R code for the chapter, in section order. |

## Code structure

### PSA simulations and EVPI

Simulates 10,000 PSA draws for a two-treatment decision model with population size 9,685. Computes:
- `maxnb[s] = max(nb0[s], nb1[s])` — the net benefit under perfect information for each simulation.
- `ol[s] = maxnb[s] - nb1[s]` — the Opportunity Loss: the amount lost by not having perfect information, given that t=2 is the overall most cost-effective option.
- EVPI via two equivalent routes: `E[maxnb] - max(E[nb0], E[nb1])` and `E[OL]`, confirming they are numerically equal.

A summary table (`tab2`) adds the VoI columns (NB*, OL) to the standard PSA table, styled with `tinytable` so the row-wise optimal NB is italicised and the Average row appears as a footer.

### EVPI curve (Vaccine dataset, BCEA)

Loads the `Vaccine` dataset, constructs the `bcea()` object and plots the EVPI vs willingness-to-pay curve using `BCEA::evi.plot()`. The "kink" at `m$kstar` is the break-even point where the optimal treatment switches. The small EVPI (≤ £2.50 per person) indicates limited decision uncertainty in this dataset.

### EVPPI (regression-based, BCEA + voi)

`BCEA::createInputs()` extracts the full S×Q PSA simulation matrix from the Vaccine model object. `BCEA::evppi()` uses a regression approximation (GAM by default) to estimate the EVPPI for three focal parameters. The EVPI–EVPPI overlay shows the EVPI as an upper bound. `info.rank()` produces the single-parameter EVPPI ranking for all parameters — a guide to research priorities, with the caveat that VoI is non-additive.

### BART illustration

Demonstrates how BART recovers a non-linear (sinusoidal) relationship between a single "covariate" x and a noisy outcome y, as an analogy for the function g(phi) = E_psi|phi[NB_t(theta)] in the EVPPI regression. `dbarts::bart()` fits 200 trees by default. The single-tree piecewise constant fit is visualised alongside the smooth ensemble prediction. BART is one of the methods available in the `voi` package for EVPPI computation, particularly effective when the number of focal parameters is large.

A demonstration table (`tab3`) shows the EVPPI regression framing: NB_t(theta) as the outcome `y` and PSA draws for phi as the covariates `X`, with nuisance parameters psi shown separately.

### EVSI (Chemotherapy example, voi package)

Uses the pre-computed EVSI object `chemo_evsi_or` bundled with `voi` (from the Chemotherapy example in Heath et al. 2024). The proposed study is a two-arm trial with binary outcome (side-effects), updating only the log OR between arms. Two plots are produced:
- EVSI vs WTP threshold for a range of sample sizes, with EVPI and EVPPI overlaid to show the hierarchy of bounds.
- EVSI vs sample size for selected WTP values, illustrating diminishing marginal returns.

### Power-based vs ENBS-based sample size

**Standard power calculation.** `power.prop.test()` computes the power curve for a two-proportion test over a grid of sample sizes (n = 0, 10, ..., 1500 per arm), assuming theta1 = 0.47 (background side-effect rate) and OR = 0.55 (chemotherapy effect). The 0.8 power threshold gives approximately 190 individuals per arm.

**ENBS analysis.** `voi::enbs()` takes the EVSI object and adds study cost assumptions (setup £5M–£10M, per-patient £28K–£42K, population 46,000, time horizon 10 years) to compute the Expected Net Benefit of Sampling at each sample size. The resulting curve shows that:
- The power-based n ≈ 190 per arm is economically suboptimal.
- The ENBS is maximised at a larger sample size (roughly 450 per arm).
- ENBS provides credible intervals (50% and 95%) for the net benefit estimate.

The key insight: standard power-based designs are framed in terms of statistical error rates; ENBS-based designs are framed in terms of monetary value, providing a more coherent basis for HTA decisions.

## Dependencies

### R packages

```r
library(tidyverse)   # data manipulation and ggplot2
library(BCEA)        # bcea(), evi.plot(), evppi(), info.rank(), createInputs(),
                     #   Vaccine dataset
library(voi)         # evpi(), enbs(), chemo_evsi_or, chemo_cea_501 datasets
library(dbarts)      # bart() for the BART illustration
library(tinytable)   # tt(), group_tt(), style_tt() for the PSA tables
library(scales)      # dollar_format() for the ENBS y-axis
```

`voi` is available from CRAN:

```r
install.packages("voi")
```

The package website with worked examples and vignettes:
[https://chjackson.github.io/voi/](https://chjackson.github.io/voi/)

The companion textbook for VoI methods:

> Heath, A., Kunst, N., & Jackson, C. (2024). *The Value of Information for Healthcare Decision Making*. CRC Press.

### LaTeX / tikz note

The VoI process cycle diagram (fig-voi-process) and the EVSI scheme (fig-evsi-scheme, imported as a PNG) have no corresponding R code in this script.
