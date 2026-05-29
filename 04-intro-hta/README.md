# Chapter: Introduction to Health Technology Assessment

## Overview

This chapter introduces the statistical and decision-theoretic framework for health technology assessment (HTA). It defines the core building blocks of a cost-effectiveness analysis — the statistical model, economic model, uncertainty analysis and decision analysis — and develops the main quantities used throughout the rest of the book: QALYs, the ICER, the monetary net benefit, the expected incremental benefit (EIB) and probabilistic sensitivity analysis (PSA). It also contrasts the frequentist ("two-stage") and Bayesian ("integrated") approaches to uncertainty propagation, and introduces the CEAC, CEAF and CEEF as graphical summaries of decision uncertainty.

## Scripts

| File | Description |
|---|---|
| `intro-hta.R` | All R code for the chapter, in section order. |

## Code structure

**QALY computation.** Demonstrates the trapezoid-rule calculation of QALYs from a sequence of utility scores measured at irregular time points. `ggbrace::stat_brace()` is used to annotate the interval width (δ_j) and the trapezoidal height ((u_j + u_{j-1})/2) on the plot.

**ICER and cost-effectiveness plane.** Manually constructs the cost-effectiveness plane for four pairwise comparisons of a reference intervention against B, C, D and E, shading the NE and SW quadrants to highlight where the ICER alone is insufficient to determine dominance.

**PSA simulation scaffolding.** Generates 10,000 illustrative PSA iterations for two interventions using `rnorm()` as stand-ins for posterior simulations. Computes the expected net benefit for each intervention, the EIB and constructs the summary table (tbl-psa) using `tinytable`, with italic and bold formatting to indicate the row-wise and overall optimal intervention.

**Cost-effectiveness plane with BCEA.** Loads the built-in `Vaccine` dataset from the `BCEA` package and constructs four versions of the CE plane: (i) the base-case ICER as a single point; (ii) the full PSA cloud; (iii) the sustainability area at k=25,000; and (iv) the sustainability area at k=10,000. The proportion of PSA points inside the sustainability area gives an intuitive measure of decision confidence.

**CEAC — two interventions.** Uses `BCEA::ceac.plot()` on the `Vaccine` dataset to show the Cost-Effectiveness Acceptability Curve: Pr(IB(θ) > 0) as a function of willingness to pay k.

**Multi-intervention analysis — Smoking dataset.** Loads the built-in `Smoking` dataset (four interventions, t=4 as reference) and produces:
- Pairwise CEACs via `ceac.plot()`.
- Per-intervention probability of being most cost-effective, computed by `multi.ce()` and plotted manually with `directlabels::geom_dl()` for inline labelling.
- Cost-Effectiveness Acceptability Frontier (CEAF) via `ceaf.plot()`: the upper envelope of the per-intervention curves, tracing the optimal intervention at each k.
- Cost-Effectiveness Efficiency Frontier (CEEF) via `ceef.plot()`: the Pareto frontier in the average cost–effectiveness space.

## Dependencies

### R packages

```r
library(tidyverse)      # data manipulation and ggplot2
library(BCEA)           # bcea(), ceac.plot(), ceaf.plot(), ceef.plot(),
                        #   multi.ce(), ceplane.plot(); Vaccine and Smoking data
library(ggbrace)        # stat_brace() for annotating the QALY diagram
library(directlabels)   # geom_dl() for inline curve labels
library(tinytable)      # tt(), style_tt(), format_tt() for the PSA table
library(scales)         # comma formatter for x-axis labels
```

`ggbrace` may need to be installed from GitHub:

```r
# install.packages("remotes")
remotes::install_github("NicolasH2/ggbrace")
```

### Data

All datasets used (`Vaccine`, `Smoking`) are bundled with the `BCEA` package and loaded via `data()`. No external files are required.

### LaTeX / tikz note

Three process-diagram figures in this chapter (the HTA decision process, the two-stage frequentist diagram and the Bayesian integrated diagram) are rendered from `.tex` source files via the `tikz` engine. They have no corresponding R code and are not reproduced in this script.
