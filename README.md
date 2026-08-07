# Bayesian modelling in health technology assessment — R code examples

This repository contains the R code accompanying the book

> **Bayesian modelling in health technology assessment**  
> Gianluca Baio  
> [gianluca.statistica.it/books/online/bmhta/](https://gianluca.statistica.it/books/online/bmhta/)

Each folder corresponds to a chapter and contains an annotated `.R` script and a `README.md` describing the code in detail. The scripts are extracted and reorganised from the book's source files, with comments added to explain the statistical intent of each code block.

---

## Repository structure

| Folder | Chapter |
|---|---|
| [`01-bayesian-reasoning`](01-bayesian-reasoning/) | Introduction to Bayesian reasoning |
| [`02-bayesian-computation`](02-bayesian-computation/) | Learning from data — Bayesian computation |
| [`03-bayesian-software`](03-bayesian-software/) | Bayesian software |
| [`04-intro-hta`](04-intro-hta/) | Introduction to health technology assessment |
| [`05-ild`](05-ild/) | Cost-effectiveness analysis with individual-level data |
| [`06-ald`](06-ald/) | Aggregated level data and evidence synthesis |
| [`07-nma`](07-nma/) | Indirect treatment comparisons (network meta-analysis) |
| [`08-survival`](08-survival/) | Survival analysis in HTA |
| [`09-markov-models`](09-markov-models/) | Markov models |
| [`10-missing-data`](10-missing-data/) | Missing data and structural values in HTA |
| [`11-itc`](11-itc/) | Population adjustment |
| [`12-voi`](12-voi/) | Value of information |

A `data/` folder at the root contains datasets shared across chapters.

---

## How to use this repository

Clone the repository and open `bmhta.Rproj` in RStudio. The project root is set as the working directory, so all relative paths (e.g. `data/ild/10TT_synth_280921.csv`) will resolve correctly without modification.

```r
# Install any missing packages before running a chapter script, e.g.
install.packages(c("tidyverse", "R2jags", "BCEA", "survHE", "survextrap",
                   "voi", "missingHE", "multinma"))
```

Each chapter folder has its own `README.md` with the full list of required packages, any external software (JAGS, Stan) and data dependencies specific to that chapter.

---

## Common dependencies

All chapters require R (≥ 4.2) and the `tidyverse`. The table below lists the main additional packages by chapter area.

| Area | Packages |
|---|---|
| Core Bayesian computation | `R2jags`, `rstan` (via `survHEhmc`) |
| Health economic evaluation | `BCEA` |
| Companion utilities | `bmhe` (see below) |
| Survival modelling | `survHE`, `survHEhmc`, `survextrap` |
| Evidence synthesis | `R2jags`, `R2OpenBUGS` |
| Population adjustment | `multinma`, `outstandR` |
| Missing data | `missingHE` |
| Value of information | `voi`, `BCEA` |

### The `bmhe` package

`bmhe` is the companion R package for this book. It provides utility functions used throughout (`stats()`, `coefplot()`, `logit()`, `ilogit()`, `gammaPar()`, `lognPar()`, etc.). Install it from the r-universe:

```r
install.packages(
  "bmhe",
  repos = c("https://giabaio.r-universe.dev", "https://cloud.r-project.org")
)
```

### External software

Several chapters require **JAGS** (chapters 2–10) and/or **Stan** (chapters 8, 9, 12):

- JAGS: [https://mcmc-jags.sourceforge.io](https://mcmc-jags.sourceforge.io)
- Stan (via rstan): [https://mc-stan.org/rstan/](https://mc-stan.org/rstan/)

---

## A note on figures

Many figures in the book use `dev: "tikz"` to render axis labels and annotations as LaTeX maths. The `.R` scripts in this repository use the default ggplot2 device, so axis labels involving maths (e.g. `$\\theta$`) will appear as plain strings rather than typeset maths when run interactively. The underlying computations are identical.

---

## Chapter summaries

**01 — Introduction to Bayesian reasoning.** Frequentist significance testing vs Bayesian induction; Bayes' theorem applied to Covid-19 testing; prior elicitation for the drug example (Beta and logit-Normal); Monte Carlo simulation and convergence; prior predictive (forward) sampling.

**02 — Bayesian computation.** Prior specification and vague priors; conjugate Beta-Binomial analysis; Covid-19 vaccine efficacy (Pfizer/BioNTech); penalised complexity priors; Jeffreys' priors; Gibbs sampling for a semi-conjugated Normal model; MCMC diagnostics (R-hat, ESS, MCSE, autocorrelation); Hamiltonian Monte Carlo illustration.

**03 — Bayesian software.** Running JAGS from R via `R2jags`; post-processing MCMC output; default differences between BUGS and JAGS; the zero-trick for non-standard distributions.

**04 — Introduction to HTA.** QALY computation; ICER and cost-effectiveness plane; PSA simulation; CEAC, CEAF and CEEF using `BCEA`; multi-intervention analysis with the Smoking dataset.

**05 — Individual-level data.** End-to-end Bayesian CEA of the 10TT trial; QALY computation with discounting; Normal/Normal independent and MCF models; Gamma/Gamma MCF; posterior predictive g-computation; DIC, WAIC and LOO-CV; structural PSA via model averaging with `BCEA`.

**06 — Aggregated level data and evidence synthesis.** Prior checks for logit-scale priors; no-pooling, complete-pooling and partial-pooling models for the magnesium meta-analysis; DIC comparison and model averaging; shrinkage visualisation; multiparameter evidence synthesis for influenza prophylaxis; cost-effectiveness analysis.

**07 — Network meta-analysis.** Evidence network diagram; data wrangling to long format; fixed-effect and random-effect NMA in JAGS; direct vs indirect evidence consistency checks; heterogeneity assessment; Half-Cauchy priors; DIC comparison; cost-effectiveness analysis using life-years gained and `BCEA`.

**08 — Survival analysis.** Kaplan-Meier estimation; Bayesian parametric survival models (Weibull, Gompertz, log-Normal, Generalised F) via `survHE`/Stan; hazard and cumulative hazard plots; PSA via posterior survival bands; DIC model selection; separate modelling by arm; M-spline models via `survextrap`; external aggregated data for extrapolation anchoring.

**09 — Markov models.** Dirichlet simplex visualisations; four-state HIV Markov model (Multinomial-Dirichlet, relative risk evidence synthesis, Weibull treatment waning, matrix-algebra simulation, discounted costs); three-state cancer model linked to survival analysis via `survHE`; Markov trace and LYG.

**10 — Missing data.** MCAR, MAR and MNAR mechanisms; simulation study of bias; multiple imputation illustrated on synthetic data; the MenSS trial with Normal-Normal, Beta-Gamma (rescaled) and Beta-Gamma hurdle models under MAR.

**11 — Population adjustment.** Effect modifiers and collapsibility (RD, RR, OR); population-average conditional vs marginal estimands; Monte Carlo vs Quasi-Monte Carlo integration (Halton sequence); overview of MAIC, STC, ML-NMR (`multinma`), g-computation and MIM (`outstandR`).

**12 — Value of information.** EVPI from PSA simulations (Opportunity Loss); EVPI curve via `BCEA`; EVPPI via regression approximation (GAM, GP, BART); Info Rank plot; EVSI and ENBS using the `voi` package; standard power-based vs ENBS-based sample size for the Chemotherapy example.

---

## Citation

If you use this code in your own work, please cite the book:

> Baio, G. (2026). *Bayesian modelling in health technology assessment*. Chapman & Hall / CRC Press.

The book is available at the [CRC website](https://www.routledge.com/Bayesian-Models-in-Health-Technology-Assessment/Baio/p/book/9781041297383). More information can be found at the [book webpage](https://gianluca.statistica.it/books/bmhta/).

---

## License

Code in this repository is released under the [MIT License](https://opensource.org/licenses/MIT). Please see the book for the full theoretical exposition accompanying each script.
