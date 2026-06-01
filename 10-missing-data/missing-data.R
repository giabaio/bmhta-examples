#' ---
#' title: "Missing data and structural values in HTA"
#' desc:  "Bayesian and frequentist approaches to missing data in HTA.
#'         Covers: complete case analysis and the implicit MCAR assumption;
#'         a simulation study comparing MNAR/MAR/MCAR bias; multiple imputation
#'         illustrated on simulated data; JAGS templates for MAR and MNAR
#'         models; and the MenSS trial example with three joint models for
#'         bivariate (e,c) outcomes under MAR: Normal-Normal MCF, Beta-Gamma
#'         MCF (rescaled), and Beta-Gamma hurdle for structural ones."
#' ---

library(tidyverse)

# ==============================================================================
# SECTION: Complete case analysis — implicit MCAR assumption
# ==============================================================================

#' Standard R regression silently performs complete case analysis (CCA),
#' which is only valid under MCAR.  The code below demonstrates this:
#' 100 observations are simulated, 20% are deleted at random, and
#' lm() automatically drops the missing rows.

x = rnorm(100, 2, 5)
y = .5 + 1.5 * x + rnorm(100)

#' Remove 20 values at random (MCAR: missingness unrelated to y or x)
y[sample(1:100, 20)] = NA

#' lm() performs CCA implicitly — note the "20 observations deleted" message.
#' With MCAR and large enough n, estimates for beta0=0.5 and beta1=1.5
#' should be approximately unbiased.
summary(lm(y ~ x))


# ==============================================================================
# SECTION: Simulation study — bias under MNAR, MAR and MCAR
# ==============================================================================

#' Runs S=1000 simulation repetitions for each missingness mechanism.
#' True model: y = 0.5 + 1.5*x + Normal(0,1), n=2000.
#'
#' In each repetition:
#'   (a) Generate full data (y, x).
#'   (b) Apply the MoM to create missing values.
#'   (c) Run CCA (lm) and save the coefficient estimates.
#'
#' MNAR: pi = ilogit(-0.5 + 0.35*x + 0.83*y) — depends on y itself.
#'   The large coefficient on y (0.83) creates substantial bias in CCA.
#' MAR:  pi = ilogit(0.17 + 0.35*x + 0*y)   — depends only on observed x.
#'   delta0 tweaked so overall missingness rate stays comparable to MNAR.
#' MCAR: pi = ilogit(0.5  + 0*x    + 0*y)   — constant probability.

# --- MNAR simulation ---
mnar = list(); pmis = numeric(); S = 1000; set.seed(140873)
for (s in 1:S) {
  x    = rnorm(2000, 2, 5)
  y    = 0.5 + 1.5 * x + rnorm(2000)
  pi   = bmhe::ilogit(-0.5 + 0.35 * x + 0.83 * y)   # MNAR: y influences pi
  m    = rbinom(2000, 1, pi)
  y[m == 1] = NA
  mnar[[s]] = (lm(y ~ x) |> summary())$coefficients[, "Estimate"]
  pmis[s]   = sum(m == 1) / 2000
}
#' Key output: beta0 and beta1 should deviate from (0.5, 1.5) due to MNAR bias.
mnar |> bind_rows() |> bmhe::stats()

# --- MAR simulation ---
mar = list(); set.seed(140873)
for (s in 1:S) {
  x  = rnorm(2000, 2, 5)
  y  = 0.5 + 1.5 * x + rnorm(2000)
  pi = bmhe::ilogit(0.17 + 0.35 * x + 0 * y)   # MAR: y has no direct effect
  m  = rbinom(2000, 1, pi)
  y[m == 1] = NA
  mar[[s]] = (lm(y ~ x) |> summary())$coefficients[, "Estimate"]
}
#' Key output: estimates should be approximately unbiased (beta0≈0.5, beta1≈1.5).
ttt = mar |> bind_rows() |> bmhe::stats()
cat("MAR analysis\n"); ttt

# --- MCAR simulation ---
mcar = list(); set.seed(140873)
for (s in 1:S) {
  x  = rnorm(2000, 2, 5)
  y  = 0.5 + 1.5 * x + rnorm(2000)
  pi = bmhe::ilogit(0.5 + 0 * x + 0 * y)   # MCAR: completely random
  m  = rbinom(2000, 1, pi)
  y[m == 1] = NA
  mcar[[s]] = (lm(y ~ x) |> summary())$coefficients[, "Estimate"]
}
#' Key output: estimates approximately unbiased, with slightly higher precision than MAR.
ttt = mcar |> bind_rows() |> bmhe::stats()
cat("MCAR analysis\n"); ttt

# --- No-missing baseline ---
nomis = list(); set.seed(140873)
for (s in 1:S) {
  x  = rnorm(2000, 2, 5)
  y  = 0.5 + 1.5 * x + rnorm(2000)
  pi = rep(0, 2000)   # no missingness
  m  = rbinom(2000, 1, pi)
  y[m == 1] = NA
  nomis[[s]] = (lm(y ~ x) |> summary())$coefficients[, "Estimate"]
}
#' Key output: tightest estimates — reference for the precision loss from missing data.
ttt = nomis |> bind_rows() |> bmhe::stats()
cat("No-missing analysis\n"); ttt


# ==============================================================================
# SECTION: Multiple imputation — illustrative example
# ==============================================================================

#' Demonstrates the key mechanics of MI:
#'   1. Generate a "population" of N=1000 with 5 covariates.
#'   2. Draw a sample, artificially delete 20% under MCAR.
#'   3. Fit a predictive model to the observed data.
#'   4. Draw R=5 imputations from the predictive distribution.
#'   5. Visualise the observed data, the fitted regression and the R imputations.

library(arm)   # for sigma.hat()

# --- Simulate population ---
N       = 1000
x0.pop  = rep(1, N)
x1.pop  = rnorm(N)
x2.pop  = rbinom(N, 1, .6)
x3.pop  = rpois(N, 3)
x4.pop  = round(rnorm(N, 40, 8))
center  = function(x) x - mean(x)
X.pop   = cbind(x0.pop, x1.pop, x2.pop, x3.pop, x4.pop)

beta      = c(.86, -1.3, .94, 2.4, 1.8)
sigma     = 6
linpred   = X.pop %*% beta
y.pop     = rnorm(N, linpred, sigma)
mu        = mean(y.pop)
sigma     = sqrt((N - 1) * sd(y.pop) / N)

# --- Draw sample with planned dropout inflation ---
w     = 1.5   # maximum 95% CI half-width
S_ss  = 4     # assumed SD for sample size calculation
y.bar = 0
nstar = 1 / ((w / 4 / S_ss)^2 + (1 / N))  # optimal sample size
p     = .2                                   # 20% dropout
n.d   = nstar / (1 - p)                     # inflated sample size

y    = sample(y.pop,    n.d, FALSE)
x0   = sample(X.pop[,1], n.d, FALSE)
x1   = sample(X.pop[,2], n.d, FALSE)
x2   = sample(X.pop[,3], n.d, FALSE)
x3   = sample(X.pop[,4], n.d, FALSE)
x4   = sample(X.pop[,5], n.d, FALSE)
X    = cbind(x0, x1, x2, x3, x4)

#' MCAR: drop out at rate p (independent of y or x)
m         = rbinom(n.d, 1, p)
y.true    = y
y[m == 1] = NA
data.mcar = data.frame(y, X, m, y.true)

# --- Fit predictive model to observed data ---
mod        = lm(y ~ x4, data = data.mcar)
mu_hat     = cbind(x0, x4) %*% mod$coefficients
sigma.hat  = sigma.hat(mod)

#' Generate R=5 imputations: for each missing y[i], draw from Normal(mu_hat[i], sigma_hat)
R = 5
y.sim = list()
for (k in 1:R) {
  y.sim[[k]] = y
  for (i in 1:length(y)) {
    y.sim[[k]][i] = ifelse(
      is.na(y.sim[[k]][i]),
      rnorm(1, mu_hat[i], sigma.hat),
      y.sim[[k]][i]
    )
  }
}

#' Combine into a tidy tibble: observed + 5 imputed datasets
ysim = tibble(y = y.sim[[1]], sim = 1) |>
  bind_rows(tibble(y = y.sim[[2]], sim = 2)) |>
  bind_rows(tibble(y = y.sim[[3]], sim = 3)) |>
  bind_rows(tibble(y = y.sim[[4]], sim = 4)) |>
  bind_rows(tibble(y = y.sim[[5]], sim = 5)) |>
  mutate(x = rep(x4, 5))

dat = tibble(x = x4, y = y, sim = "obs") |>
  bind_rows(tibble(x = x4, y = y, sim = "1") |> dplyr::filter(is.na(y))) |>
  bind_rows(tibble(x = x4, y = y, sim = "2") |> dplyr::filter(is.na(y))) |>
  bind_rows(tibble(x = x4, y = y, sim = "3") |> dplyr::filter(is.na(y))) |>
  bind_rows(tibble(x = x4, y = y, sim = "4") |> dplyr::filter(is.na(y))) |>
  bind_rows(tibble(x = x4, y = y, sim = "5") |> dplyr::filter(is.na(y))) |>
  rowwise() |>
  mutate(
    y = case_when(
      sim != "obs" ~ rnorm(1, cbind(1, x) %*% mod$coefficients, sigma.hat),
      TRUE ~ y
    )
  ) |>
  mutate(type = case_when(sim == "obs" ~ 16, TRUE ~ 18))

# --- Panel (a): partially observed data ---
dat |> dplyr::filter(sim == "obs") |>
  ggplot(aes(x, y)) + geom_point() +
  ylab("(Partially) observed response") + xlab("Observed covariate") +
  annotate("text", x = (dat |> dplyr::filter(is.na(y)))$x, -Inf,
           label = "?", vjust = -1) +
  ylim(range(dat$y, na.rm = TRUE))

# --- Panel (b): fitted linear regression ---
dat |> dplyr::filter(sim == "obs") |>
  ggplot(aes(x, y)) + geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ylab("Observed response") + xlab("Observed covariate") +
  ylim(range(dat$y, na.rm = TRUE))

# --- Panel (c): R=5 imputations ---
dat |>
  ggplot(aes(x, y, col = as.factor(sim), shape = as.factor(sim))) +
  geom_point(size = 1.6) +
  scale_color_manual(
    name   = "",
    values = c(rep("#649eff", 5), "black"),
    labels = c(paste("Imputation", 1:5), "Observed")
  ) +
  scale_shape_manual(
    name   = "",
    values = c(21, 22, 23, 24, 25, 19),
    labels = c(1:5, "Observed")
  ) +
  guides(
    color = guide_legend(
      override.aes = list(shape = c(21, 22, 23, 24, 25, 19))
    ),
    shape = "none"
  ) +
  xlab("Observed covariate") + ylab("Observed/Imputed outcome")



# ==============================================================================
# SECTION: MenSS trial — data loading and descriptive plots
# ==============================================================================

#' The MenSS (Men's Safer Sex) pilot RCT evaluates a digital STI intervention
#' vs standard care.  QALYs and costs are measured at baseline and at 3, 6
#' and 12 months.  Missing data are substantial, especially in the active arm
#' (only 23% complete cases in the intervention arm).
#'
#' Data loaded from a pre-processed .rds file; see gabrio2019 for details.

data = readRDS("data/missing-data/data_MenSS.rds")

#' Assemble a single long-format tibble.
datalist = list(
  N    = data$N1 + data$N2,
  trt  = c(rep(1, data$N1), rep(2, data$N2)),
  eff  = c(data$eff1, data$eff2),
  q0   = c(data$q0_1, data$q0_2),
  cost = c(data$cost1, data$cost2)
) |> as_tibble()

#' Distribution of QALYs by arm.
#' The spike at exactly 1 (highlighted in orange) reflects the structural ones:
#' participants in perfect health who account for a non-trivial share of the data.
datalist |>
  mutate(treat = factor(case_when(trt == 1 ~ "Control", trt == 2 ~ "Intervention"))) |>
  ggplot(aes(eff, fill = eff == 1)) +
  geom_histogram(col = "black") + xlab("QALYs") + ylab("") + theme_bw() +
  facet_grid(. ~ treat) +
  scale_fill_manual(values = c("#1F77B4", "#FF7F0E")) +
  theme(legend.position = "none")

#' Distribution of costs by arm.
datalist |>
  mutate(treat = factor(case_when(trt == 1 ~ "Control", trt == 2 ~ "Intervention"))) |>
  ggplot(aes(cost)) +
  geom_histogram(fill = "grey", col = "black") +
  facet_grid(. ~ treat) + xlab("Costs") + ylab("")

#' Summary table of observed proportions at each time point.
mensstab = tibble(
  time    = c("Baseline", "3 months", "6 months", "12 months", "Complete cases"),
  outcome = c(
    "utilities",
    "utilities and costs", "utilities and costs",
    "utilities and costs", "utilities and costs"
  ),
  observed1 = c("72 (96%)", "34 (45%)", "35 (47%)", "43 (57%)", "27 (44%)"),
  observed2 = c("72 (86%)", "23 (27%)", "23 (27%)", "36 (43%)", "19 (23%)")
)
colnames(mensstab) = c("Time", "Type of outcome", "observed (%)", "observed (%)")
mensstab |> tinytable::tt() |> tinytable::style_tt(i = 6, j = 1, bold = TRUE)


# ==============================================================================
# SECTION: MenSS trial — three models (pre-computed results)
# ==============================================================================

#' The three JAGS models are run in a companion script
#' (data/missing-data/menss-analysis.R) and results saved to
#' MenSS_JAGS_models.rds.  See that script for the full JAGS model code.
#'
#' Briefly, all three assume MAR for effects.  The three approaches differ
#' in how they handle the QALYs distribution:
#'
#' Model 1 — Normal-Normal MCF
#'   e_i ~ Normal(phi_ei, tau_e)
#'   phi_ei = alpha0 + alpha1*Trt[i] + alpha2*(u0[i] - ubar)
#'   c_i ~ Normal(phi_ci, tau_c)
#'   phi_ci = beta0 + beta1*Trt[i] + beta2*(e[i] - mu_e[Trt[i]])
#'   Conceptually simple but unconstrained: predicted QALYs can exceed 1.
#'
#' Model 2 — Beta-Gamma MCF (rescaled)
#'   logit(phi_ei) = alpha0 + alpha1*Trt[i] + alpha2*(u0[i] - ubar)
#'   e_i ~ Beta(nu_e * phi_ei, nu_e * (1 - phi_ei))
#'   Costs modelled with log-Normal or Gamma MCF.
#'   To avoid the structural spike at 1, all QALYs are shifted by -eps=0.01.
#'   Correctly bounded but the epsilon shift introduces slight downward bias
#'   and ignores the distinct "healthy" subgroup.
#'
#' Model 3 — Beta-Gamma MCF with hurdle for structural ones
#'   For each individual: d_i ~ Bernoulli(gamma_i) indicates perfect health.
#'   Those with d_i=1 are modelled as structural ones (QALY=1 fixed).
#'   Those with d_i=0 are modelled with Beta, no rescaling needed.
#'   Population average: mu_e[t] = (1 - gamma_bar[t]) * mu_e_lt1[t] + gamma_bar[t]
#'   This correctly captures the mixture structure.

#' Load pre-computed JAGS results.
objs = readRDS("data/missing-data/MenSS_JAGS_models.rds")

#' Posterior predictive plot for all three models.
#' Dots = posterior predicted mean for missing individuals.
#' Segments = 95% predictive interval.
#' Orange crosses = observed data.
objs$cp1 |> mutate(model = "Normal-Normal") |>
  bind_rows(objs$cp2 |> mutate(model = "Beta-Gamma (rescaled)")) |>
  bind_rows(objs$cp3 |> mutate(model = "Beta-Gamma hurdle")) |>
  mutate(
    model = factor(
      model,
      levels = c("Normal-Normal", "Beta-Gamma (rescaled)", "Beta-Gamma hurdle")
    ),
    col = case_when(col == "black" ~ "#FF7F0E", TRUE ~ "#1F77B4")
  ) |>
  ggplot(aes(mean, ID, shape = factor(shp), col = col)) +
  geom_point() +
  scale_shape_manual(values = c(16, 4)) +
  theme(legend.position = "none") +
  geom_linerange(aes(xmin = `2.5%`, xmax = `97.5%`)) +
  facet_grid(rows = vars(model), cols = vars(arm), scales = "free_x") +
  xlab("QALYs") + ylab("Individuals") +
  coord_flip() + geom_vline(xintercept = 1, linetype = "dashed") +
  scale_color_identity() +
  scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0.0, 1.2))

#' Alternative: overlay only Beta-Gamma rescaled and hurdle (omit Normal-Normal).
objs$cp1 |> mutate(model = "Normal-Normal") |>
  bind_rows(objs$cp2 |> mutate(model = "Beta-Gamma\n (rescaled)")) |>
  bind_rows(objs$cp3 |> mutate(model = "Beta-Gamma\n hurdle")) |>
  mutate(
    model = factor(
      model,
      levels = c("Normal-Normal", "Beta-Gamma\n (rescaled)", "Beta-Gamma\n hurdle")
    ),
    shp = factor(shp)
  ) |>
  dplyr::filter(model != "Normal-Normal") |>
  ggplot(aes(x = num, y = mean, color = model, shape = shp)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`),
                 position = position_dodge(width = 0.5)) +
  scale_shape_manual(values = c(16, 4)) +
  scale_color_manual(values = c(
    "Beta-Gamma\n (rescaled)" = "#1F77B4",
    "Beta-Gamma\n hurdle"     = "#FF7F0E"
  )) +
  facet_wrap(. ~ arm, scales = "free_x", ncol = 1) +
  xlab("Individuals") + ylab("QALYs") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  coord_cartesian(ylim = c(0.3, 1.0)) +
  theme(
    legend.position   = "none",
    axis.ticks.x      = element_blank(),
    axis.text.x       = element_blank()
  )
