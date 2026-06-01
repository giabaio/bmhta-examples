# Chapter: Population adjustment

## Overview

This chapter addresses the problem of indirect treatment comparisons (ITCs) when the patient populations of different trials are not exchangeable. It introduces the conceptual framework of population adjustment, discusses effect modifiers and the collapsibility of common treatment effect measures, and reviews the main methods used in HTA practice: MAIC, STC, multilevel network meta-regression (ML-NMR), parametric g-computation, and multiple imputation marginalisation (MIM). The chapter is predominantly methodological; the computational code provides two self-contained illustrations: the collapsibility example and a Monte Carlo vs Quasi-Monte Carlo comparison.

## Scripts

| File | Description |
|---|---|
| `itc.R` | All R code for the chapter, in section order. |

## Code structure

The chapter's R code falls into three sections. The two main applied methods (ML-NMR and MIM) are not implemented directly here — they are available via the `multinma` and `outstandR` packages respectively, with worked examples in those packages' documentation (see Dependencies below).

### Collapsibility — contingency table

Builds the fictional Treatment vs Placebo trial from Chancellor et al. (adapted), stratified by smoking status. The table is rendered with `tinytable::group_tt()` for column grouping. Three treatment effect measures are then computed stratum-by-stratum and in aggregate:

- **Risk Difference (RD)**: stratum-specific values are equal (no effect modification) and the crude overall RD equals the stratum-weighted average — RD is *collapsible*.
- **Odds Ratio (OR)**: stratum-specific ORs are equal (no effect modification on this scale) but the crude OR differs from the log-scale weighted average — OR is *not collapsible*.
- **Risk Ratio (RR)**: stratum-specific RRs differ — smoking IS an effect modifier on the RR scale, even though it is not on the RD or OR scales.

This illustrates why the choice of effect measure matters for indirect comparisons: naive pooling of ORs or HRs across populations with different covariate distributions will produce biased estimates even when there is no confounding.

### Population-average conditional and marginal estimands

Using the same contingency table, computes the two estimands that ML-NMR and g-computation target respectively, both on the logit scale:

- **Conditional estimand** (`d_BA`): weighted sum of stratum-specific logit contrasts, where weights are the covariate distribution in the target population. This is the average individual-level treatment effect on the linear predictor scale — equivalent to the ML-NMR estimand (eq-mlnmr-cond-releff).

- **Marginal estimand** (`delta_BA`): first marginalise the predicted probabilities over the covariate distribution to get overall means `bar_mu_B` and `bar_mu_A`, then take their logit contrast. This is the estimand targeted by g-computation and MIM. The marginal means coincide numerically with the crude estimates from the "Total" column of the table.

Because the logit is non-collapsible, `d_BA != delta_BA`. The difference quantifies the "aggregation bias" that arises when conditional effects are reported in place of marginal ones for economic decision-making.

Uses `bmhe::logit()` throughout.

### Monte Carlo vs Quasi-Monte Carlo integration

Illustrates the computational motivation for QMC within ML-NMR. The target integral is the area of a circle of radius 0.3 centred at (0.5, 0.5), approximated using n=100 function evaluations:

- **Monte Carlo**: n points drawn from Uniform(0,1) independently. Pseudo-random sampling creates visible clustering and gaps, leading to a less accurate estimate.
- **Quasi-Monte Carlo**: n points from a Halton sequence (base 2 for x, base 3 for y), implemented from scratch as `halton(n, base)`. The Halton sequence is a low-discrepancy deterministic sequence that fills the unit square more evenly, generally achieving O(n^{-1}) convergence for smooth integrands vs O(n^{-1/2}) for MC.

Two ggplot panels show the sampling points overlaid on the circle region, with triangles for points where f=1 and dots for f=0. The plot title reports the analytic value, the estimate and the absolute error for each method.

## Dependencies

### R packages

```r
library(tidyverse)   # data manipulation and ggplot2
library(bmhe)        # logit(), ilogit()
library(tinytable)   # tt(), group_tt(), style_tt() for the contingency table
```

`bmhe` is the companion package for this book:

```r
install.packages(
  "bmhe",
  repos = c("https://giabaio.r-universe.dev", "https://cloud.r-project.org")
)
```

### Packages for the methods described (not implemented in this script)

**ML-NMR**: the `multinma` package (Phillippo 2020+), available from CRAN. Worked examples at [https://dmphillippo.github.io/multinma/](https://dmphillippo.github.io/multinma/).

```r
install.packages("multinma")
```

**G-computation and MIM**: the `outstandR` package (Remiro-Azócar 2024+), available from GitHub. Worked examples at [https://statisticshealtheconomics.github.io/outstandR/](https://statisticshealtheconomics.github.io/outstandR/).

```r
remotes::install_github("StatisticsHealthEconomics/outstandR")
```

Both packages use Stan (via `rstan` or `cmdstanr`) under the hood.

### LaTeX / tikz note

The network structure diagram, the g-computation imputation-stage diagram (fig-mim1), the MIM analysis-stage diagram (fig-mim2) and the missingness mechanism DAGs are rendered from `.tex` source files via the tikz engine and have no corresponding R code.
