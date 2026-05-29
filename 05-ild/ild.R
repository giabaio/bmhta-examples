#' ---
#' title: "Cost-effectiveness analysis with individual-level data"
#' desc:  "End-to-end Bayesian cost-effectiveness analysis of the 10TT trial
#'         (synthetic data). Covers: QALY computation with discounting; three
#'         JAGS models for joint (e,c) outcomes (Normal/Normal independent,
#'         Normal/Normal MCF, Gamma/Gamma MCF); posterior predictive g-computation;
#'         model selection via DIC/pD/pV; WAIC and LOO-CV; and cost-effectiveness
#'         analysis using BCEA including structural PSA via model averaging."
#' ---

library(tidyverse)
library(R2jags)
library(BCEA)


# ==============================================================================
# SECTION: Load data
# ==============================================================================

#' Loads the synthetic 10TT dataset and keeps only complete cases.
#' The original trial randomised 537 obese adults to a brief weight-loss
#' intervention (10 Top Tips leaflet) vs standard care; this synthetic version
#' mimics the key structure.  After dropping missing records, n=167.
ttt = read.csv("data/ild/10TT_synth_280921.csv") |>
  as_tibble() |>
  drop_na()


# ==============================================================================
# SECTION: QALY computation with discounting
# ==============================================================================

#' Computes discounted QALYs for each individual using the trapezoid rule
#' applied to HRQL utility scores measured at baseline (0), 3, 6, 12, 18
#' and 24 months.
#'
#' Discounting reduces the value of utilities that accrue in the future
#' (NICE recommends d=3.5% per year).  Here, the baseline month (0) is
#' assigned to year 1 (pmax ensures no division by (1+d)^0 = 1 at baseline).

disc = function(x, year, disc_rate = 0.035) {
  #' Applies a discount factor to an outcome measured at a given year.
  #' x:         the outcome value to discount
  #' year:      the year at which the measurement occurs (baseline = year 1)
  #' disc_rate: annual discount rate (default 3.5%)
  x / ((1 + disc_rate)^(year - 1))
}

#' Pivots the wide-format utility columns to long format (one row per
#' individual-month), applies the discount and then computes the trapezoid
#' AUC for each individual.  The result is the sum of all trapezoids divided
#' by 12 to convert from months to years.
df_qol = ttt |>
  select(id, arm, contains("qol")) |>
  pivot_longer(
    contains("qol"),
    names_to        = "month",
    names_prefix    = "qol_",
    names_transform = list(month = as.integer),
    values_to       = "qol"
  ) |>
  # Assign each measurement to a calendar year; baseline stays in year 1
  mutate(year = pmax(1, ceiling(month / 12))) |>
  mutate(qol_disc = disc(qol, year)) |>
  group_by(id) |>
  mutate(
    # Duration of the interval between two consecutive measurements (months)
    delta = month - dplyr::lag(month, default = 0),
    # Sum of utilities at the start and end of the interval (trapezoid height)
    du    = qol_disc + dplyr::lag(qol_disc, default = 0),
    # Area of this trapezoid (half base × height)
    auc   = du * delta / 2
  ) |>
  # Sum trapezoids and convert from months to years
  summarise(qaly = sum(auc) / 12) |>
  ungroup()

#' Merges computed QALYs back into the main dataset and relabels the arm
#' variable for readability.
ttt = ttt |>
  left_join(df_qol, by = c("id" = "id")) |>
  mutate(
    Trt = arm + 1,
    arm = factor(arm),
    arm = case_when(arm == "0" ~ "Control", TRUE ~ "Treatment")
  )

#' Visual check: distributions of discounted QALYs and total costs by arm.
#' QALYs are expected to be left-skewed (most individuals maintain good
#' health over 2 years); costs are typically right-skewed.
ttt |>
  ggplot(aes(qaly)) +
  geom_histogram(fill = "grey", col = "black") +
  facet_grid(. ~ as.factor(arm)) +
  xlab("QALYs") + ylab("")

ttt |>
  ggplot(aes(totalcost)) +
  geom_histogram(fill = "grey", col = "black") +
  facet_grid(. ~ as.factor(arm)) +
  xlab("Total costs (GBP)") + ylab("")


# ==============================================================================
# SECTION: Prior visualisation
# ==============================================================================

#' Plots the three families of prior distributions used across the models.
#'
#' Regression coefficients: Normal(0, 100) — widely spread, effectively vague.
#' sigma_e (effects SD):    Exponential(5.75) — PC prior with Pr(sigma_e>0.8)≈0.01.
#' sigma_c (costs SD):      Exponential(0.35) — PC prior with Pr(sigma_c>2)≈0.5,
#'                          on the £1000 scale (costs are divided by 1000 before
#'                          fitting to keep priors on a sensible scale).

ggplot() +
  stat_function(fun = dnorm, args = list(mean = 0, sd = 100)) +
  xlim(-500, 500) + xlab("") + ylab("")

ggplot() +
  stat_function(fun = dexp, args = list(rate = 5.75)) +
  xlim(0, 2) + xlab("") + ylab("")

ggplot() +
  stat_function(fun = dexp, args = list(rate = 0.35)) +
  xlim(0, 10) + xlab("£ x 1000") + ylab("")


# ==============================================================================
# SECTION: Data preparation for JAGS
# ==============================================================================

#' Prepares the data objects passed to all three JAGS models.
#' Costs are scaled to £1000 to keep their range compatible with the
#' Normal(0, 100) priors on regression coefficients.
#' u0star is the centred baseline utility (zero-mean); centering aids
#' interpretation of the intercept and helps convergence.
e      = ttt$qaly
c      = ttt$totalcost / 1000        # rescale to £1000 units
Trt    = ttt$Trt
u0star = scale(ttt$qol_0, scale = FALSE) |> as.numeric()
N      = nrow(ttt)

data   = list(e = e, c = c, Trt = Trt, u0star = u0star, N = N)

#' We use a single named list `model` to store all JAGS output objects,
#' making it easy to compare models and feed them into BCEA later.
model = list()


# ==============================================================================
# SECTION: Model 1 — Normal/Normal independent
# ==============================================================================

#' Models effects and costs as independent Normal linear regressions.
#' This is the "standard" approach and serves as a baseline for comparison.
#'
#' Model for effects (QALYs):
#'   e_i ~ Normal(phi_ei, tau_e[Trt[i]])
#'   phi_ei = alpha0 + alpha1*(Trt[i]-1) + alpha2*u0star[i]
#'
#' Model for costs:
#'   c_i ~ Normal(phi_ci, tau_c[Trt[i]])
#'   phi_ci = beta0 + beta1*(Trt[i]-1)
#'
#' Population averages (t=1 control, t=2 treatment):
#'   mu.e[t] = alpha0 + alpha1*(t-1)
#'   mu.c[t] = 1000*(beta0 + beta1*(t-1))   [rescale back to GBP]
#'
#' Priors:
#'   Regression coefficients: Normal(0, 0.0001) [= Normal(mean=0, precision=0.0001)]
#'   sigma.e[t], sigma.c[t]:  Exponential PC priors (arm-specific)

nn_indep = function() {
  for (i in 1:N) {
    e[i] ~ dnorm(phi.e[i], tau.e[Trt[i]])
    phi.e[i] <- alpha0 + alpha1 * (Trt[i] - 1) + alpha2 * u0star[i]
    c[i] ~ dnorm(phi.c[i], tau.c[Trt[i]])
    phi.c[i] <- beta0 + beta1 * (Trt[i] - 1)
  }
  for (t in 1:2) {
    mu.e[t] <- alpha0 + alpha1 * (t - 1)
    mu.c[t] <- 1000 * (beta0 + beta1 * (t - 1))
  }
  alpha0 ~ dnorm(0, 0.0001)
  alpha1 ~ dnorm(0, 0.0001)
  alpha2 ~ dnorm(0, 0.0001)
  beta0  ~ dnorm(0, 0.0001)
  beta1  ~ dnorm(0, 0.0001)
  for (t in 1:2) {
    sigma.e[t] ~ dexp(5.75)      # PC prior: Pr(sigma_e > 0.8) ≈ 0.01
    tau.e[t]   <- pow(sigma.e[t], -2)
    sigma.c[t] ~ dexp(0.35)      # PC prior: Pr(sigma_c > 2) ≈ 0.5
    tau.c[t]   <- pow(sigma.c[t], -2)
  }
}

model$nn_indep = jags(
  data               = data,
  parameters.to.save = c(
    "mu.e", "mu.c", "alpha0", "alpha1", "alpha2",
    "beta0", "beta1", "sigma.e", "sigma.c"
  ),
  inits      = NULL,
  n.chains   = 2, n.iter = 5000, n.burnin = 3000, n.thin = 1,
  DIC        = TRUE,
  model.file = nn_indep
)

print(model$nn_indep, digits = 3, interval = c(0.025, 0.5, 0.975))


# ==============================================================================
# SECTION: Model 2 — Normal/Normal MCF (marginal/conditional factorisation)
# ==============================================================================

#' Extends the independence model by adding a regression of costs on effects
#' within the conditional model for costs.  This allows the model to capture
#' correlation between e and c without a bivariate Normal specification.
#'
#' The key addition is the coefficient beta2 in the cost linear predictor:
#'   phi_ci = beta0 + beta1*(Trt[i]-1) + beta2*(e[i] - mu.e[Trt[i]])
#'
#' Because both e and c are Normal, this is algebraically equivalent to the
#' Seemingly Unrelated Regression (SUR) / bivariate Normal model, and the
#' correlation and marginal cost SD can be derived analytically:
#'   sigma.c[t]^2 = lambda.c[t]^2 + sigma.e[t]^2 * beta2^2
#'   rho[t]       = beta2 * sigma.e[t] / sigma.c[t]
#'
#' lambda.c[t] is the *conditional* SD for costs given effects; sigma.c[t]
#' is the resulting *marginal* SD for costs.

nn_mcf = function() {
  for (i in 1:N) {
    # Marginal model for effects
    e[i] ~ dnorm(phi.e[i], tau.e[Trt[i]])
    phi.e[i] <- alpha0 + alpha1 * (Trt[i] - 1) + alpha2 * u0star[i]
    # Conditional model for costs given effects
    c[i] ~ dnorm(phi.c[i], tau.c[Trt[i]])
    phi.c[i] <- beta0 + beta1 * (Trt[i] - 1) + beta2 * (e[i] - mu.e[Trt[i]])
  }
  for (t in 1:2) {
    mu.e[t] <- alpha0 + alpha1 * (t - 1)
    mu.c[t] <- 1000 * (beta0 + beta1 * (t - 1))
  }
  alpha0 ~ dnorm(0, 0.0001)
  alpha1 ~ dnorm(0, 0.0001)
  alpha2 ~ dnorm(0, 0.0001)
  beta0  ~ dnorm(0, 0.0001)
  beta1  ~ dnorm(0, 0.0001)
  beta2  ~ dnorm(0, 0.0001)
  for (t in 1:2) {
    sigma.e[t]  ~ dexp(5.75)
    tau.e[t]    <- pow(sigma.e[t], -2)
    lambda.c[t] ~ dexp(0.35)      # conditional SD for costs
    tau.c[t]    <- pow(lambda.c[t], -2)
    # Marginal SD and correlation derived analytically
    rho[t]     <- beta2 * sigma.e[t] / sigma.c[t]
    sigma.c[t] <- sqrt(pow(lambda.c[t], 2) + pow(sigma.e[t], 2) * pow(beta2, 2))
  }
}

model$nn_mcf = jags(
  data               = data,
  parameters.to.save = c(
    "mu.e", "mu.c", "alpha0", "alpha1", "alpha2",
    "beta0", "beta1", "beta2", "sigma.e", "sigma.c", "lambda.c", "rho"
  ),
  inits      = NULL,
  n.chains   = 2, n.iter = 5000, n.burnin = 3000, n.thin = 1,
  DIC        = TRUE,
  model.file = nn_mcf
)

print(model$nn_mcf, digits = 3, interval = c(0.025, 0.5, 0.975))

#' Alternative: computing rho and sigma.c post-hoc in R rather than inside JAGS.
#' This is equivalent to the JAGS-computed values and can be more efficient
#' for complex models where minimising the number of monitored nodes matters.
lambda.c = model$nn_mcf$BUGSoutput$sims.list$lambda.c
sigma.e  = model$nn_mcf$BUGSoutput$sims.list$sigma.e
beta1    = model$nn_mcf$BUGSoutput$sims.list$beta1

sigma.c = rho = matrix(NA, nrow = nrow(lambda.c), ncol = ncol(lambda.c))
for (t in 1:2) {
  sigma.c[, t] = sqrt(lambda.c[, t]^2 + (sigma.e[, t]^2 * beta1^2))
  rho[, t]     = beta1 * sigma.e[, t] / sigma.c[, t]
}
colnames(sigma.c) = c("sigma.c[1]", "sigma.c[2]")
colnames(rho)     = c("rho[1]", "rho[2]")
cbind(sigma.c, rho) |> bmhe::stats() |> round(digits = 3)

#' Coefficient plot comparing the two Normal/Normal models.
#' The two models produce very similar estimates; this is expected when the
#' correlation is small, because the MCF model reduces to the independence model.
toplot1 = bmhe::coefplot(model$nn_indep, parameter = c("alpha", "mu.e", "sigma.e"))$data
toplot2 = bmhe::coefplot(model$nn_mcf,   parameter = c("alpha", "mu.e", "sigma.e", "rho"))$data

toplot1 |> ggplot(aes(mean, Parameter)) +
  geom_linerange(aes(xmin = low, xmax = upp), position = position_nudge(y = 0.1)) +
  geom_point(position = position_nudge(y = 0.1)) + theme_bw() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Interval estimate") +
  geom_linerange(data = toplot2, aes(xmin = low, xmax = upp),
                 position = position_nudge(y = -.1), col = "red") +
  geom_point(data = toplot2, position = position_nudge(y = -0.1), col = "red")

toplot1 = bmhe::coefplot(model$nn_indep, parameter = c("beta", "mu.c", "sigma.c"))$data
toplot1[3, 2:5] = toplot1[3, 2:5] / 1000
toplot1[4, 2:5] = toplot1[4, 2:5] / 1000
toplot2 = bmhe::coefplot(model$nn_mcf, parameter = c("beta", "mu.c", "sigma.c"))$data
toplot2[4, 2:5] = toplot2[4, 2:5] / 1000
toplot2[5, 2:5] = toplot2[5, 2:5] / 1000

toplot1 |> ggplot(aes(mean, Parameter)) +
  geom_linerange(aes(xmin = low, xmax = upp), position = position_nudge(y = 0.1)) +
  geom_point(position = position_nudge(y = 0.1)) + theme_bw() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Interval estimate") +
  geom_linerange(data = toplot2, aes(xmin = low, xmax = upp),
                 position = position_nudge(y = -.1), col = "red") +
  geom_point(data = toplot2, position = position_nudge(y = -0.1), col = "red")


# ==============================================================================
# SECTION: Model 3 — Gamma/Gamma MCF
# ==============================================================================

#' Replaces the Normal sampling distribution with Gamma for both effects and
#' costs, appropriate because:
#'   (a) total costs are strictly positive and right-skewed;
#'   (b) QALYs over a 2-year horizon are bounded above and left-skewed.
#'
#' Data transformation for QALYs: estar[i] = 3 - e[i].
#'   This flips the left skew into a right skew, making the Gamma distribution
#'   an appropriate model.  Because the transformation is linear,
#'   E[e] = 3 - E[estar], so population average QALYs can be recovered as
#'   mu.e[t] = 3 - mustar.e[t].
#'
#' The Gamma distribution is parameterised by shape (nu) and rate (gamma=nu/phi),
#' where phi = E[Y] is the mean.  Log-linear predictors are used to ensure
#' the mean is positive:
#'   log(phi.e[i]) = alpha0 + alpha1*(Trt[i]-1) + alpha2*u0star[i]
#'   log(phi.c[i]) = beta0 + beta1*(Trt[i]-1) + beta2*(estar[i]-mustar.e[Trt[i]])
#'
#' PC prior on shape parameters: Exponential(0.15), encoding Pr(nu>30)≈0.01.
#' Costs are still scaled to £1000; mu.c is rescaled back at the end.

#' Add the rescaled QALY variable to the data list.
data$estar = 3 - data$e

gg_mcf = function() {
  for (i in 1:N) {
    # Marginal model for rescaled effects
    estar[i] ~ dgamma(nu.e[Trt[i]], gamma.e[Trt[i], i])
    gamma.e[Trt[i], i] <- nu.e[Trt[i]] / phi.e[i]
    log(phi.e[i]) <- alpha0 + alpha1 * (Trt[i] - 1) + alpha2 * u0star[i]
    # Conditional model for costs given rescaled effects
    c[i] ~ dgamma(nu.c[Trt[i]], gamma.c[Trt[i], i])
    gamma.c[Trt[i], i] <- nu.c[Trt[i]] / phi.c[i]
    log(phi.c[i]) <- beta0 + beta1 * (Trt[i] - 1) +
      beta2 * (estar[i] - mustar.e[Trt[i]])
  }
  for (t in 1:2) {
    mustar.e[t] <- exp(alpha0 + alpha1 * (t - 1))  # mean on rescaled scale
    mu.e[t]     <- 3 - mustar.e[t]                 # back-transform to original
    mu.c[t]     <- 1000 * exp(beta0 + beta1 * (t - 1))  # rescale to GBP
  }
  alpha0 ~ dnorm(0, 0.0001)
  alpha1 ~ dnorm(0, 0.0001)
  alpha2 ~ dnorm(0, 0.0001)
  beta0  ~ dnorm(0, 0.0001)
  beta1  ~ dnorm(0, 0.0001)
  beta2  ~ dnorm(0, 0.0001)
  for (t in 1:2) {
    nu.e[t] ~ dexp(0.15)   # PC prior: Pr(nu.e > 30) ≈ 0.01
    nu.c[t] ~ dexp(0.15)   # PC prior: Pr(nu.c > 30) ≈ 0.01
  }
}

model$gg_mcf = jags(
  data               = data,
  parameters.to.save = c(
    "mu.e", "mustar.e", "mu.c",
    "alpha0", "alpha1", "alpha2",
    "beta0",  "beta1",  "beta2",
    "nu.e",   "nu.c"
  ),
  inits      = NULL,
  n.chains   = 2, n.iter = 5000, n.burnin = 3000, n.thin = 1,
  DIC        = TRUE,
  model.file = gg_mcf
)

print(model$gg_mcf, digits = 3, interval = c(0.025, 0.5, 0.975))

#' Comparison of Gamma and log-Normal densities at matching mean (8) and SD (5),
#' illustrating the key difference: the log-Normal has heavier tails.
#' This matters practically: Gamma is typically the safer choice for costs because
#' it assigns less prior probability to very extreme values.
set.seed(14)
x = seq(0, 100, .1)
m = 8; s = 5
g  = dgamma(x, shape = m^2 / s^2, rate = m / s^2)
ln = dlnorm(x, bmhe::lognPar(m, s)$mulog, bmhe::lognPar(m, s)$sigmalog)

toplot = tibble(
  x     = x, d = g,
  model = paste0(
    "Gamma$(\\nu,\\gamma)$: mean=", m, ", sd=", s, ", median=",
    format(median(rgamma(10000, shape = m^2 / s^2, rate = m / s^2)), digits = 3, nsmall = 3)
  )
) |>
  bind_rows(
    tibble(
      x     = x, d = ln,
      model = paste0(
        "log-Normal$(\\eta,\\lambda)$: mean=", m, ", sd=", s, ", median=",
        format(
          median(rlnorm(10000, bmhe::lognPar(m, s)$mulog, bmhe::lognPar(m, s)$sigmalog)),
          digits = 3, nsmall = 3
        )
      )
    )
  )

# Panel (a): full range
tibble(x = 12, y = 0.065) |>
  ggplot(aes(x, y)) + xlab("$y$") + ylab("") +
  geom_line(data = toplot, aes(x, d, col = model, linetype = model)) +
  theme(
    legend.position        = c(0.65, 0.85),
    legend.background      = element_rect(fill = "transparent"),
    legend.title           = element_blank()
  ) + xlim(0, 40) +
  scale_color_manual(values    = c("#000000", "#1F77B4"), name = "") +
  scale_linetype_manual(values = c("solid",   "dashed"),  name = "")

# Panel (b): right tail — log-Normal has noticeably more mass
tibble(x = 20, y = 0.065) |>
  ggplot(aes(x, y)) + xlab("$y$") + ylab("") +
  geom_line(data = toplot, aes(x, d, col = model, linetype = model)) +
  theme(
    legend.position        = c(0.65, 0.85),
    legend.background      = element_rect(fill = "transparent"),
    legend.title           = element_blank()
  ) + xlim(50, 100) + ylim(0, .0001) +
  scale_color_manual(values    = c("#000000", "#1F77B4"), name = "") +
  scale_linetype_manual(values = c("solid",   "dashed"),  name = "")


# ==============================================================================
# SECTION: Posterior predictive g-computation
# ==============================================================================

#' When the transformation from the model scale to the natural scale is
#' non-linear, we cannot simply apply the inverse function to posterior means.
#' Instead we use posterior predictive Monte Carlo:
#'
#'   For each MCMC draw s and each arm t:
#'     1. Simulate N values of estar from Gamma(nu.e[s,t], nu.e[s,t]/mustar.e[s,t])
#'     2. Back-transform each draw: e = 3 - estar
#'     3. Compute the mean across the N draws -> mu[s, t]
#'
#' The resulting matrix mu is a sample from the posterior distribution of the
#' population average QALYs on the original scale.
#' For linear transformations (as here, e = 3 - estar), this is equivalent
#' to the direct calculation mu.e[t] = 3 - mustar.e[t] already in the JAGS code.
#' The general procedure is essential for non-linear transformations.

nu.e    = model$gg_mcf$BUGSoutput$sims.list$nu.e
mustar.e = model$gg_mcf$BUGSoutput$sims.list$mustar.e
rate    = nu.e / mustar.e

N_mc = 4000   # MC draws for the predictive step
mu   = matrix(NA, nrow = nrow(rate), ncol = 2)

for (i in 1:nrow(nu.e)) {
  for (t in 1:2) {
    estar    = rgamma(N_mc, shape = nu.e[i, t], rate = rate[i, t])
    mu[i, t] = mean(3 - estar)
  }
}

colnames(mu) = c("mu.e[1]", "mu.e[2]")
mu |> bmhe::stats()
# Results should closely match model$gg_mcf$BUGSoutput$summary rows "mu.e[1]" and "mu.e[2]"


# ==============================================================================
# SECTION: Manual computation of pD and DIC
# ==============================================================================

#' Demonstrates how to compute p_D (the BUGS-style DIC penalty) manually from
#' R, for the Gamma/Gamma MCF model.
#'
#' p_D = Dbar - Dhat, where:
#'   Dbar = E[D(theta)] = posterior mean deviance (stored by JAGS)
#'   Dhat = D(E[theta]) = deviance evaluated at the posterior mean parameters
#'
#' This is useful because JAGS by default computes pV = Var[D]/2, which is
#' invariant to reparameterisation but numerically different from pD.
#' To get pD directly from R2jags, use pD=TRUE in the jags() call (see below).

# --- Effects component of the likelihood ---

#' Extracts posterior mean parameters for the effects model.
nu.e     = c(
  model$gg_mcf$BUGSoutput$summary["nu.e[1]", "mean"],
  model$gg_mcf$BUGSoutput$summary["nu.e[2]", "mean"]
)
alpha0 = model$gg_mcf$BUGSoutput$summary["alpha0", "mean"]
alpha1 = model$gg_mcf$BUGSoutput$summary["alpha1", "mean"]
alpha2 = model$gg_mcf$BUGSoutput$summary["alpha2", "mean"]

#' Reconstructs the linear predictor and rate parameter for each individual,
#' then evaluates the Gamma log-likelihood at the posterior mean parameters.
trt        = data$Trt - 1
u0         = data$u0star
log.phi.e  = alpha0 + alpha1 * trt + alpha2 * u0

mustar.e   = numeric()
for (t in 1:2) {
  mustar.e[t] = exp(alpha0 + alpha1 * (t - 1))
}
gamma.e = numeric()
for (i in 1:data$N) {
  gamma.e[i] = nu.e[data$Trt[i]] / exp(log.phi.e[i])
}
# Individual deviance contributions from the effects
lik.e = -2 * dgamma(data$estar, shape = nu.e[data$Trt], rate = gamma.e, log = TRUE)

# --- Costs component of the likelihood ---

nu.c  = c(
  model$gg_mcf$BUGSoutput$summary["nu.c[1]", "mean"],
  model$gg_mcf$BUGSoutput$summary["nu.c[2]", "mean"]
)
beta0 = model$gg_mcf$BUGSoutput$summary["beta0", "mean"]
beta1 = model$gg_mcf$BUGSoutput$summary["beta1", "mean"]
beta2 = model$gg_mcf$BUGSoutput$summary["beta2", "mean"]

#' Note the regression of costs on (estar - mustar.e) within each arm,
#' which is the MCF adjustment for correlation.
log.phi.c = beta0 + beta1 * trt + beta2 * (data$estar - mustar.e[data$Trt]) 
gamma.c   = numeric()
for (i in 1:data$N) {
  gamma.c[i] = nu.c[data$Trt[i]] / exp(log.phi.c[i])
}
lik.c = -2 * dgamma(data$c, shape = nu.c[data$Trt], rate = gamma.c, log = TRUE)

# --- Assemble pD and DIC ---

dbar = model$gg_mcf$BUGSoutput$summary["deviance", "mean"]   # Dbar
dhat = sum(lik.e) + sum(lik.c)                               # Dhat
pD   = dbar - dhat
DIC  = dbar + pD   # equivalent to dhat + 2*pD
cat("pD =", round(pD, 2), "\nDIC =", round(DIC, 1), "\n")


# ==============================================================================
# SECTION: Model selection — re-run all models with pD=TRUE
# ==============================================================================

#' Store existing pV and DIC values from the initial runs before re-running.
pv  = c(
  model$nn_indep$BUGSoutput$pV,
  model$nn_mcf$BUGSoutput$pV,
  model$gg_mcf$BUGSoutput$pV
)
dic = c(
  model$nn_indep$BUGSoutput$DIC,
  model$nn_mcf$BUGSoutput$DIC,
  model$gg_mcf$BUGSoutput$DIC
)

#' Adding pD=TRUE causes R2jags to call rjags::dic.samples() after the main
#' MCMC run to compute the BUGS-style penalty pD.  The result is stored in
#' model$...$BUGSoutput$pD and model$...$BUGSoutput$DIC2.
model$nn_indep = jags(
  data               = data,
  parameters.to.save = c(
    "mu.e", "mu.c", "alpha0", "alpha1", "alpha2", "beta0", "beta1",
    "sigma.e", "sigma.c"
  ),
  inits      = NULL,
  n.chains   = 2, n.iter = 5000, n.burnin = 3000, n.thin = 1,
  DIC        = TRUE,
  model.file = nn_indep,
  pD         = TRUE
)
print(model$nn_indep, digits = 3, interval = c(0.025, 0.5, 0.975))

model$nn_mcf = jags(
  data               = data,
  parameters.to.save = c(
    "mu.e", "mu.c", "alpha0", "alpha1", "alpha2",
    "beta0", "beta1", "beta2", "sigma.e", "sigma.c", "lambda.c", "rho"
  ),
  inits      = NULL,
  n.chains   = 2, n.iter = 5000, n.burnin = 3000, n.thin = 1,
  DIC        = TRUE,
  model.file = nn_mcf,
  pD         = TRUE
)

model$gg_mcf = jags(
  data               = data,
  parameters.to.save = c(
    "mu.e", "mustar.e", "mu.c",
    "alpha0", "alpha1", "alpha2", "beta0", "beta1", "beta2",
    "nu.e", "nu.c"
  ),
  inits      = NULL,
  n.chains   = 2, n.iter = 5000, n.burnin = 3000, n.thin = 1,
  DIC        = TRUE,
  model.file = gg_mcf,
  pD         = TRUE
)

#' Summary table of DIC-based model selection statistics.
#' Rule of thumb: ΔDIC < 2 = equivalent fit; 3–7 = less support;
#' > 10 = negligible support for the worse-fitting model.
tab = tibble(
  Model = c("Normal/Normal independent", "Normal/Normal MCF", "Gamma/Gamma MCF"),
  pv    = pv |> round(2),
  DIC   = dic |> round(),
  pd    = c(
    model$nn_indep$BUGSoutput$pD |> round(2),
    model$nn_mcf$BUGSoutput$pD  |> round(2),
    model$gg_mcf$BUGSoutput$pD  |> round(2)
  ),
  DIC2  = c(
    model$nn_indep$BUGSoutput$DIC2 |> round(),
    model$nn_mcf$BUGSoutput$DIC2  |> round(),
    model$gg_mcf$BUGSoutput$DIC2  |> round()
  )
)
colnames(tab) = c(
  "Model", "pV", "DIC (pV)", "pD", "DIC (pD)"
)
tab |>
  tinytable::tt() |>
  tinytable::style_tt(j = 2:5, align = "c")


# ==============================================================================
# SECTION: WAIC and LOO-CV (Gamma/Gamma MCF)
# ==============================================================================

#' WAIC and LOO-CV are alternatives to DIC based on the posterior predictive
#' distribution rather than the deviance.  They require the individual
#' log-likelihood contributions to be monitored.
#'
#' The JAGS model below is identical to gg_mcf() but adds:
#'   log.lik[i] = logdensity.gamma(estar[i], ...) + logdensity.gamma(c[i], ...)
#' where logdensity.gamma() is a built-in JAGS function for the Gamma log-density.

gg_mcf_waic = function() {
  for (i in 1:N) {
    estar[i] ~ dgamma(nu.e[Trt[i]], gamma.e[Trt[i], i])
    gamma.e[Trt[i], i] <- nu.e[Trt[i]] / phi.e[i]
    log(phi.e[i]) <- alpha0 + alpha1 * (Trt[i] - 1) + alpha2 * u0star[i]
    c[i] ~ dgamma(nu.c[Trt[i]], gamma.c[Trt[i], i])
    gamma.c[Trt[i], i] <- nu.c[Trt[i]] / phi.c[i]
    log(phi.c[i]) <- beta0 + beta1 * (Trt[i] - 1) +
      beta2 * (estar[i] - mustar.e[Trt[i]])
    # Individual log-likelihood = sum of effects and costs contributions
    log.lik[i] <- logdensity.gamma(estar[i], nu.e[Trt[i]], gamma.e[Trt[i], i]) +
                  logdensity.gamma(c[i],     nu.c[Trt[i]], gamma.c[Trt[i], i])
  }
  for (t in 1:2) {
    mustar.e[t] <- exp(alpha0 + alpha1 * (t - 1))
    mu.e[t]     <- 3 - mustar.e[t]
    mu.c[t]     <- exp(beta0 + beta1 * (t - 1))
  }
  alpha0 ~ dnorm(0, 0.0001)
  alpha1 ~ dnorm(0, 0.0001)
  alpha2 ~ dnorm(0, 0.0001)
  beta0  ~ dnorm(0, 0.0001)
  beta1  ~ dnorm(0, 0.0001)
  beta2  ~ dnorm(0, 0.0001)
  for (t in 1:2) {
    nu.e[t] ~ dexp(0.15)
    nu.c[t] ~ dexp(0.15)
  }
}

m_waic = jags(
  data               = data,
  parameters.to.save = c("log.lik"),
  inits              = NULL,
  n.chains           = 2, n.iter = 5000, n.burnin = 3000, n.thin = 1,
  DIC                = TRUE,
  model.file         = gg_mcf_waic
)

#' loo::waic() takes a matrix of log-likelihood values (n.sims × n_obs) and
#' returns the WAIC and its penalty pW.
loo::waic(m_waic$BUGSoutput$sims.list$log.lik)

#' LOO-CV is more robust than WAIC when some individual pW values are large
#' (a warning from waic() suggests this).  loo() uses Pareto smoothed
#' importance sampling to stabilise the leave-one-out estimates.
loo::loo(m_waic$BUGSoutput$sims.list$log.lik)


# ==============================================================================
# SECTION: DIC weights for model averaging
# ==============================================================================

#' When several models have non-trivial DIC differences, we can combine their
#' predictions as a weighted average rather than selecting a single winner.
#' The weight for model h decays steeply as its DIC moves away from the minimum:
#'   w_h = exp(-0.5 * ΔDIC_h) / sum_h exp(-0.5 * ΔDIC_h)
#'
#' Plotted here as a function of ΔDIC for the case of two models (one being
#' the best-fitting).  The weight effectively vanishes at ΔDIC ≈ 10, consistent
#' with the rule-of-thumb for DIC comparisons.

w = function(x) {
  for (i in 1:length(x)) {
    w = (exp(-.5 * c(x, 0)) / sum(exp(-.5 * c(x, 0))))[1]
  }
  w
}
tibble(x = seq(0, 100, .1)) |>
  rowwise() |>
  mutate(y = w(x)) |>
  ggplot(aes(x, y)) + geom_line() + xlim(0, 10) +
  xlab("Absolute difference |min DIC - DIC_h|") +
  ylab("Weight associated with model h")


# ==============================================================================
# SECTION: Cost-effectiveness analysis using BCEA
# ==============================================================================

#' Passes the MCMC posterior simulations for population average effects and
#' costs directly to BCEA.  Each bcea() call constructs the cost-effectiveness
#' plane, computes the ICER, EIB, CEAC etc. and stores them in an object that
#' all BCEA plotting functions can consume.

interventions = c("Standard", "Active intervention")
ref           = 2    # intervention arm is the reference for comparison

# --- Normal/Normal independent ---
eff  = model$nn_indep$BUGSoutput$sims.list$mu.e
cost = model$nn_indep$BUGSoutput$sims.list$mu.c
m_nn_indep = bcea(eff = eff, cost = cost, interventions = interventions, ref = ref)

# --- Normal/Normal MCF ---
eff  = model$nn_mcf$BUGSoutput$sims.list$mu.e
cost = model$nn_mcf$BUGSoutput$sims.list$mu.c
m_nn_mcf = bcea(eff = eff, cost = cost, interventions = interventions, ref = ref)

# --- Gamma/Gamma MCF ---
eff  = model$gg_mcf$BUGSoutput$sims.list$mu.e
cost = model$gg_mcf$BUGSoutput$sims.list$mu.c
m_gg_mcf = bcea(eff = eff, cost = cost, interventions = interventions, ref = ref)

#' Cost-effectiveness plane with joint posterior contour for each model.
contour2(m_nn_indep, graph = "gg")
contour2(m_nn_mcf,   graph = "gg")
contour2(m_gg_mcf,   graph = "gg")

#' CEACs for all three models, overlaid on one plot.
#' The graph="gg" option returns a ggplot object whose $data element can be
#' extracted and combined for customised multi-model plots.
p1 = ceac.plot(m_nn_indep, graph = "gg")
p2 = ceac.plot(m_nn_mcf,   graph = "gg")
p3 = ceac.plot(m_gg_mcf,   graph = "gg")

p1$data |> mutate(model = "Normal/Normal indep") |>
  bind_rows(p2$data |> mutate(model = "Normal/Normal MCF")) |>
  bind_rows(p3$data |> mutate(model = "Gamma/Gamma MCF"))   |>
  ggplot(aes(k, ceac, color = model)) + geom_line() +
  xlab("Willingness to pay") + ylim(0, 1) +
  ylab("Probability of cost effectiveness") +
  labs(title = "Cost Effectiveness Acceptability Curves", color = "Model") +
  theme(
    legend.position.inside = c(.6, .75), legend.position = "inside",
    legend.background      = element_blank()
  ) +
  scale_color_manual(values = c("#000000", "#1F77B4", "#FF7F0E"))

#' Numerical summary for the best-fitting model (Gamma/Gamma MCF).
summary(m_gg_mcf)


# ==============================================================================
# SECTION: Structural PSA via model averaging (struct.psa)
# ==============================================================================

#' Computes a DIC-weighted average of the three models' economic outputs.
#' The weights are derived from equation (eq-dic-weights) and handled internally
#' by BCEA::struct.psa().  The result is a bcea object that can be used with
#' all the usual BCEA methods.
#'
#' In this example the Gamma/Gamma MCF dominates with weight ≈ 1, so the
#' model average is nearly identical to that model alone.  Model averaging
#' is most informative when two or more models are closely matched in DIC.

effects = list(
  model$nn_indep$BUGSoutput$sims.list$mu.e,
  model$nn_mcf$BUGSoutput$sims.list$mu.e,
  model$gg_mcf$BUGSoutput$sims.list$mu.e
)
costs = list(
  model$nn_indep$BUGSoutput$sims.list$mu.c,
  model$nn_mcf$BUGSoutput$sims.list$mu.c,
  model$gg_mcf$BUGSoutput$sims.list$mu.c
)
models = list(
  "nn-indep" = model$nn_indep,
  "nn-mcf"   = model$nn_mcf,
  "gg-mcf"   = model$gg_mcf
)

m_avg = struct.psa(models, effects, costs, ref = ref, interventions = interventions)

#' Inspect the weights assigned to each model.
m_avg$w

#' Summary of the model-averaged decision.
summary(m_avg)
