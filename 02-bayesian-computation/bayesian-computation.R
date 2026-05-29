#' ---
#' title:  "Bayesian Computation"
#' desc:   "R code extracted from Chapter: Learning from data -- Bayesian computation.
#'          Covers prior specification, conjugate analysis, a worked Covid vaccine
#'          example, penalised complexity priors, Jeffreys' priors, Gibbs sampling
#'          for a semi-conjugated Normal model and MCMC diagnostics."
#' ---


# ==============================================================================
# SETUP
# ==============================================================================

#' Sets the number of Monte Carlo simulations used throughout the chapter and
#' draws a prior sample for theta and its implied predictive data y.
#' The Beta(9.2, 13.8) prior is the "drug" example running example:
#' it encodes the belief that the drug's success probability is around 0.4
#' with relatively modest spread.
nsim = 10000
alpha = 9.2
beta  = 13.8
theta = rbeta(n = nsim, shape1 = alpha, shape2 = beta)
y     = rbinom(n = nsim, size = 20, prob = theta)


# ==============================================================================
# SECTION: Vague / Uniform priors
# ==============================================================================

#' Compares the Beta(9.2, 13.8) prior against a Uniform(0,1) ("vague") prior.
#' The two panels in fig-priors show how much more spread the Uniform is --
#' it effectively says every success probability is equally plausible a priori.
#' eval=FALSE in the source because the actual figure is built in the chunks
#' below; this version is shown to readers as illustrative code.
priors = tibble(
  theta = seq(from = 0, to = 1, by = .001),
  p1    = dbeta(seq(0, 1, .001), 9.2, 13.8),
  p2    = dunif(seq(0, 1, .001), 0, 1)
)
priors |> ggplot(aes(theta, p1)) + geom_line()
priors |> ggplot(aes(theta, p2)) + geom_line()

#' Draws forward-sampling (prior predictive) simulations under both priors.
#' Comparing y1 and y2 illustrates how a vague prior implies a much wider
#' range of plausible trial outcomes -- the 95% interval for y2 covers
#' essentially the entire [0, 20] range.
theta1 = rbeta(n = nsim, shape1 = 9.2, shape2 = 13.8)
theta2 = runif(n = nsim, min = 0, max = 1)
y1 = rbinom(n = nsim, size = 20, prob = theta1)
y2 = rbinom(n = nsim, size = 20, prob = theta2)

#' Summarises the two predictive distributions side by side.
#' bmhe::stats() is a convenience wrapper from the book's companion package.
sims = cbind(y1, y2)
bmhe::stats(sims)

#' Computes the Beta-Binomial prior predictive probability mass function
#' p(y) = integral p(y|theta) p(theta) dtheta analytically, for both priors.
#' This is the marginal (prior predictive) distribution; panel heights differ
#' between the two priors because the Uniform spreads mass more evenly.
bb = function(y, n, a, b) {
  choose(n, y) * beta(y + a, n + b - y) / beta(a, b)
}
p1 = p2 = numeric()
for (i in 1:21) {
  p1[i] = bb((i - 1), 20, 9.2, 13.8)
  p2[i] = bb((i - 1), 20,   1,   1)
}
tibble(y = seq(0, 20), p = p1, label = "Beta prior") |>
  ggplot(aes(x = y, y = p)) + geom_col(fill = "grey", col = "black") +
  ylab("") + xlab("$y$")
tibble(y = seq(0, 20), p = p2, label = "Uniform prior") |>
  ggplot(aes(x = y, y = p)) + geom_col(fill = "grey", col = "black") +
  ylim(0, max(p1)) + ylab("") + xlab("$y$")

#' Posterior under the Uniform prior after observing y=15 out of n=20.
#' A Uniform prior is conjugate to the Binomial and yields Beta(y+1, n-y+1).
#' The darker horizontal segment is the Bayesian 95% credible interval;
#' the lighter one is the MLE-based 95% confidence interval for comparison.
tibble(theta = seq(0, 1, .001)) |>
  mutate(post = dbeta(theta, 15 + 1, 20 - 15 + 1)) |>
  ggplot(aes(theta, post)) + geom_line() + xlab("$\\theta$") + ylab("") +
  geom_segment(
    aes(
      x    = qbeta(0.025, 15 + 1, 20 - 15 + 1),
      xend = qbeta(0.975, 15 + 1, 20 - 15 + 1),
      y = 0, yend = 0
    ),
    col = "grey10", size = 1.1
  ) +
  geom_segment(
    aes(
      x    = binom.test(15, 20)$conf.int[1],
      xend = binom.test(15, 20)$conf.int[2],
      y = 0.05, yend = 0.05
    ),
    col = "grey80", size = 1.1
  )

#' Monte Carlo sample from the Uniform-prior posterior and a numerical summary.
post = rbeta(10000, 15 + 1, 20 - 15 + 1)
bmhe::stats(post)


# ==============================================================================
# SECTION: Eliciting priors via pseudo-data (thought experiment)
# ==============================================================================

#' Illustrates how to calibrate a prior by imagining a "pretend" dataset.
#' Starting from a Uniform(0,1) and updating with (y0, n0) pseudo-observations
#' produces a Beta(y0+1, n0-y0+1) that encodes a particular level of belief.
#' The four panels in fig-many-betas show what happens for different pseudo-data
#' choices; (y0=0, n0=0) recovers complete ignorance (= Uniform again).
y = c(8, 15, 150, 0)
n = c(20, 20, 200, 0)

tibble(theta = seq(0, 1, .001)) |>
  mutate(p1 = dbeta(theta, y[1] + 1, n[1] - y[1] + 1)) |>
  ggplot(aes(theta, p1)) + geom_line() + ylab("") + xlab("$\\theta$")

tibble(theta = seq(0, 1, .001)) |>
  mutate(p2 = dbeta(theta, y[2] + 1, n[2] - y[2] + 1)) |>
  ggplot(aes(theta, p2)) + geom_line() + ylab("") + xlab("$\\theta$") +
  geom_segment(
    aes(
      x    = qbeta(.025, y[2] + 1, n[2] - y[2] + 1),
      xend = qbeta(.975, y[2] + 1, n[2] - y[2] + 1),
      y = 0, yend = 0
    ),
    col = "blue", linewidth = 1.1
  ) +
  geom_segment(
    aes(
      x    = (y[2] + 1) / (n[2] + 2),
      xend = (y[2] + 1) / (n[2] + 2),
      y = -Inf, yend = Inf
    ),
    col = "blue", linetype = 2
  )

tibble(theta = seq(0, 1, .001)) |>
  mutate(p3 = dbeta(theta, y[3] + 1, n[3] - y[3] + 1)) |>
  ggplot(aes(theta, p3)) + geom_line() + ylab("") + xlab("$\\theta$")

tibble(theta = seq(0, 1, .001)) |>
  mutate(p4 = dbeta(theta, y[4] + 1, n[4] - y[4] + 1)) |>
  ggplot(aes(theta, p4)) + geom_line() + ylab("") + xlab("$\\theta$")


# ==============================================================================
# SECTION: Vague Normal prior (for unbounded parameters)
# ==============================================================================

#' A Normal prior with very large variance is proper (integrates to 1) yet
#' extremely diffuse, making it a sensible "vague" choice for unbounded
#' parameters such as regression coefficients.
#' The two panels simply show the same density at two different y-axis scales.
ggplot() + stat_function(fun = dnorm, args = list(mean = 0, sd = 100000)) +
  xlim(-500000, 500000) + ylim(0, .00015) + xlab("$\\theta$") + ylab("") +
  annotate("text", x =  Inf, y = 0, label = "$\\rightarrow\\infty$",   vjust = -1.5, hjust =  1.5) +
  annotate("text", x = -Inf, y = 0, label = "$-\\infty\\leftarrow$", vjust = -1.5, hjust = -.5)

ggplot() + stat_function(fun = dnorm, args = list(mean = 0, sd = 100000)) +
  xlim(-500000, 500000) + xlab("$\\theta$") + ylab("")


# ==============================================================================
# SECTION: Conjugate analysis (Beta-Binomial)
# ==============================================================================

#' Plots the prior, (rescaled) likelihood and posterior for the Beta-Binomial
#' conjugate model.  With prior Beta(a, b) and data (y, n) the posterior is
#' Beta(a+y, b+n-y) -- no MCMC required.
#' directlabels::geom_dl() places labels directly on the curves to avoid
#' a colour legend.
theta = seq(0, 1, .001)
a = 9.2
b = 13.8
r = 15
n = 20

tibble(theta = theta) |>
  mutate(d = dbeta(theta, a, b), type = "$p(\\theta)$") |>
  bind_rows(
    tibble(theta = theta, d = dbeta(theta, r + 1, n - r + 1),
           type = "$\\mathcal{L}(\\theta\\mid y=15,n=20)$")
  ) |>
  bind_rows(
    tibble(theta = theta, d = dbeta(theta, a + r, b + n - r),
           type = "$p(\\theta\\mid y=15,n=20)$")
  ) |>
  ggplot(aes(theta, d, col = type)) + geom_line() +
  directlabels::geom_dl(aes(label = type),
                        method = list("top.bumpup", hjust = 0.25, vjust = -.5)) +
  xlab("$\\theta$") + ylab("") + theme(legend.position = "none") +
  ylim(0, 5.5) +
  scale_color_manual(values = c("#000", "#000", "#000"))


# ==============================================================================
# SECTION: Covid-19 vaccine example (Pfizer-BioNTech)
# ==============================================================================

#' Replicates the Bayesian sample-size and efficacy analysis from the
#' Pfizer/BioNTech Phase II/III trial.
#'
#' The trial reformulates vaccine efficacy (VE) as theta = pi_vac/(pi_vac + pi_plac),
#' allowing a single Beta prior to be used.  The prior mean is set to correspond
#' to the minimum acceptable VE of 30%, which maps to theta = 0.4117.
#'
#' Because the two arms had slightly unequal sizes at the time of analysis,
#' the case counts are rescaled to restore the 1:1 allocation assumed by the model.
alpha.0 = 0.700102
beta.0  = 1
y = c(8, 162)
n = c(17411, 17511)

# Rescale counts to equalise arm sizes
y = y * mean(n) / n

# Conjugate update: Beta(alpha.0, beta.0) + Binomial data -> Beta(alpha.1, beta.1)
alpha.1 = alpha.0 + y[1]
beta.1  = beta.0  + y[2]

# Monte Carlo draws from the posterior for theta and the derived VE
theta = rbeta(100000, alpha.1, beta.1)
ve    = (1 - 2 * theta) / (1 - theta)

#' Visualises the posterior distributions for theta and VE.
#' The narrow histogram for VE illustrates the very high estimated efficacy.
tibble(theta) |>
  ggplot(aes(theta)) + geom_histogram(col = "black", fill = "grey") +
  xlab("$\\theta$") + ylab("")

tibble(ve) |>
  ggplot(aes(ve)) + geom_histogram(col = "black", fill = "grey") +
  xlab("VE") + ylab("")


# ==============================================================================
# SECTION: Penalised Complexity (PC) priors
# ==============================================================================

# --- PC prior for a Binomial parameter ----------------------------------------

#' Implements the PC prior for a Bernoulli probability theta.
#' The prior penalises deviation from a base-case model with theta0 = 0.5
#' (maximum entropy / equal probability), using the Kullback-Leibler divergence
#' as the distance measure.  The singularity guard (d < 1e-10) prevents
#' division by zero at theta = theta0.
pc_prior = function(theta, theta0 = 0.5, lambda = 1) {
  # theta0: baseline probability under the simpler model
  # lambda: scaling factor controlling how fast the prior decays with distance
  d     = sqrt(2 * theta * log(theta / theta0) +
               2 * (1 - theta) * log((1 - theta) / (1 - theta0)))
  deriv = ifelse(
    d < 1e-10, 1,
    abs((log(theta) - log(1 - theta) - log(theta0) + log(1 - theta0)) / d)
  )
  lambda * exp(-lambda * d) * deriv
}

#' Uses numerical integration to find the threshold u such that
#' Pr(theta > u) = alpha under the PC prior.
#' This calibrates lambda by trial-and-error: we want roughly 22% of the
#' prior mass above 0.75.
theta = seq(0, 1, .0001)
u     = seq(0.5, .9999, .01)
ints  = list()
for (i in 1:length(u)) {
  ints[[i]] = integrate(pc_prior, u[i], 1)
}
probs = numeric()
for (i in 1:length(ints)) {
  probs[i] = ints[[i]]$value
}

# More concise equivalent using lapply/unlist (as noted in the text)
ints  = lapply(u, function(x) integrate(pc_prior, x, 1))
probs = unlist(lapply(ints, `[[`, 1))

# Normalise so that probabilities sum to 1 over [0,1]
probs = probs / integrate(pc_prior, 0, 1)$value

# Find the threshold corresponding to alpha = 0.22
alpha = 0.22
u[min(which(probs <= alpha))]

#' Plots the PC prior alongside Jeffreys' Beta(0.5,0.5) and Uniform(0,1)
#' for two parameterisations, using facets.
#' Note: the plotting version of pc_prior omits the singularity guard for brevity.
pc_prior = function(theta, theta0 = 0.5, lambda = 1) {
  d     = sqrt(2 * theta * log(theta / theta0) +
               2 * (1 - theta) * log((1 - theta) / (1 - theta0)))
  deriv = abs((log(theta) - log(1 - theta) - log(theta0) + log(1 - theta0)) / d)
  lambda * exp(-lambda * d) * deriv
}

tibble(theta = theta, p = pc_prior(theta),
       type = "PC prior $\\theta_0=0.5,\\lambda=1$") |>
  bind_rows(tibble(theta = theta, p = pc_prior(theta, theta0 = 0.2, lambda = 2),
                   type = "PC prior $\\theta_0=0.2,\\lambda=2$")) |>
  bind_rows(tibble(theta = theta, p = dunif(theta, 0, 1), type = "Uniform(0,1)")) |>
  bind_rows(tibble(theta = theta, p = dbeta(theta, .5, .5), type = "Jeffreys' prior")) |>
  ggplot(aes(theta, p)) + geom_line() + ylim(0, 5) +
  xlab("$\\theta$") + ylab("Prior") +
  theme(legend.position = "none") +
  facet_wrap(. ~ type)

# --- PC prior for a standard deviation ----------------------------------------

#' Plots the type-2 Gumbel PC prior for the precision tau and the resulting
#' Exponential prior for sigma, under the constraint Pr(sigma > 1) = 0.01.
#' The Exponential prior on sigma is obtained by a change-of-variable from tau.
tau     = seq(0.01, 100, .01)
sigma_0 = 1
alpha   = 0.01
lambda  = -log(alpha) / sigma_0   # = -log(0.01)/1 ≈ 4.61

tibble(
  tau = tau,
  p   = (lambda / 2) * tau^(-3/2) * exp(-lambda * tau^(-1/2))
) |> ggplot(aes(tau, p)) + geom_line() +
  xlab(label = "$\\tau$") + ylab(label = "PC prior")

tibble(sigma = 1 / sqrt(tau)) |>
  mutate(p = dexp(sigma, lambda)) |>
  ggplot(aes(sigma, p)) + geom_line() +
  xlab("$\\sigma$") + ylab("") + xlim(0, 3)


# ==============================================================================
# SECTION: Jeffreys' prior (Binomial example)
# ==============================================================================

#' Custom legend key: displays a number label instead of the default line glyph.
#' Used in fig-jeffreys to keep the legend compact.
draw_key_numbered = function(data, params, size) {
  number_label = as.character(data$label %||% "1")
  key_rect     = draw_key_rect(data, params, size)
  key_text     = grid::textGrob(
    number_label, gp = grid::gpar(col = "black", fontsize = 10)
  )
  grid::grobTree(key_rect, key_text)
}

#' Binomial likelihood function; used in all four Jeffreys' prior panels.
lik = function(theta, y, n) {
  choose(n, y) * theta^y * (1 - theta)^(n - y)
}
theta = seq(0, 1, .001)

#' Four-panel figure comparing Jeffreys' Beta(0.5, 0.5) prior with the
#' Binomial likelihood and resulting posterior, for varying (y, n).
#' Key insight: whatever the data, the likelihood shape dominates the posterior
#' under Jeffreys' prior.

# Panel 1: y=0, n=2
y = 0; n = 2
tibble(theta = theta) |>
  mutate(d = lik(theta, y, n),
         type  = paste0("$\\mathcal{L}(\\theta\\mid y=", y, ",n=", n, ")$"),
         label = "1") |>
  bind_rows(tibble(theta = theta,
                   d     = dbeta(theta, y + 0.5, n - y + 0.5),
                   type  = paste0("Beta(", y + 0.5, ",", n - y + 0.5, ") posterior"),
                   label = "2")) |>
  bind_rows(tibble(theta = seq(0, 1, .001),
                   d     = dbeta(theta, 0.5, 0.5),
                   type  = "Jeffreys' prior",
                   label = "3")) |>
  ggplot(aes(theta, d, col = type)) +
  geom_line(key_glyph = draw_key_label) + ylim(0, 2.5) +
  guides(colour = guide_legend(override.aes = list(label = c("1", "2", "3")))) +
  directlabels::geom_dl(aes(label = label),
                        method = list("first.points", vjust = c(0, 4, 5), hjust = c(0, -5, 0))) +
  theme(legend.position.inside = c(0.55, .65), legend.position = "inside",
        legend.background = element_blank(), legend.title = element_blank()) +
  xlab("$\\theta$") + ylab("")

# Panel 2: y=2, n=2
y = 2; n = 2
tibble(theta = theta) |>
  mutate(d = lik(theta, y, n),
         type = paste0("$\\mathcal{L}(\\theta\\mid y=", y, ",n=", n, ")$"), label = "1") |>
  bind_rows(tibble(theta = theta, d = dbeta(theta, y + .5, n - y + .5),
                   type = paste0("Beta(", y + 0.5, ",", n - y + 0.5, ") posterior"), label = "2")) |>
  bind_rows(tibble(theta = seq(0, 1, .001), d = dbeta(theta, 0.5, 0.5),
                   type = "Jeffreys' prior", label = "3")) |>
  ggplot(aes(theta, d, col = type)) + geom_line(key_glyph = draw_key_label) + ylim(0, 2.5) +
  theme(legend.position = c(0.55, .65), legend.background = element_blank(),
        legend.title = element_blank()) +
  xlab("$\\theta$") + ylab("") +
  guides(colour = guide_legend(override.aes = list(label = c("1", "2", "3")))) +
  directlabels::geom_dl(aes(label = label),
                        method = list("last.points", vjust = c(0, 4, 5), hjust = c(0, 5, 0)))

# Panel 3: y=1, n=2
y = 1; n = 2
tibble(theta = theta) |>
  mutate(d = lik(theta, y, n),
         type = paste0("$\\mathcal{L}(\\theta\\mid y=", y, ",n=", n, ")$"), label = "1") |>
  bind_rows(tibble(theta = theta, d = dbeta(theta, y + .5, n - y + .5),
                   type = paste0("Beta(", y + 0.5, ",", n - y + 0.5, ") posterior"), label = "2")) |>
  bind_rows(tibble(theta = seq(0, 1, .001), d = dbeta(theta, 0.5, 0.5),
                   type = "Jeffreys' prior", label = "3")) |>
  ggplot(aes(theta, d, col = type)) + geom_line(key_glyph = draw_key_label) + ylim(0, 2.5) +
  theme(legend.position = c(0.55, .65), legend.background = element_blank(),
        legend.title = element_blank()) +
  xlab("$\\theta$") + ylab("") +
  guides(colour = guide_legend(override.aes = list(label = c("1", "2", "3")))) +
  directlabels::geom_dl(aes(label = label),
                        method = list("top.bumpup", vjust = c(-.25, -.25, 5), hjust = c(0, 0, 0)))

# Panel 4: y=2, n=10
y = 2; n = 10
tibble(theta = theta) |>
  mutate(d = lik(theta, y, n),
         type = paste0("$\\mathcal{L}(\\theta\\mid y=", y, ",n=", n, ")$"), label = "1") |>
  bind_rows(tibble(theta = theta, d = dbeta(theta, y + .5, n - y + .5),
                   type = paste0("Beta(", y + 0.5, ",", n - y + 0.5, ") posterior"), label = "2")) |>
  bind_rows(tibble(theta = seq(0, 1, .001), d = dbeta(theta, 0.5, 0.5),
                   type = "Jeffreys' prior", label = "3")) |>
  ggplot(aes(theta, d, col = type)) + geom_line(key_glyph = draw_key_label) + ylim(0, 3.5) +
  theme(legend.position = c(0.55, .65), legend.background = element_blank(),
        legend.title = element_blank()) +
  xlab("$\\theta$") + ylab("") +
  guides(colour = guide_legend(override.aes = list(label = c("1", "2", "3")))) +
  directlabels::geom_dl(aes(label = label),
                        method = list("top.bumpup", vjust = c(-.25, -.25, 5), hjust = c(0, 0, 0)))


# ==============================================================================
# SECTION: Markov Chain Monte Carlo -- Gibbs sampling
# ==============================================================================

# --- Data and prior hyper-parameters ------------------------------------------

#' Observed data for the semi-conjugated Normal model (n=30 observations).
y = c(
  1.2697,  7.7637,  2.2532,  3.4557,  4.1776,  6.4320, -3.6623,  7.7567,
  5.9032,  7.2671, -2.3447,  8.0160,  3.5013,  2.8495,  0.6467,  3.2371,
  5.8573, -3.3749,  4.1507,  4.3092, 11.7327,  2.6174,  9.4942, -2.7639,
 -1.5859,  3.6986,  2.4544, -0.3294,  0.2329,  5.2846
)

#' Prior hyper-parameters.
#' mu ~ Normal(mu_0, sigma2_0): very large variance => effectively vague prior.
#' tau ~ Gamma(alpha_0, beta_0): near-zero shape and rate => vague prior for precision.
mu_0     = 0
sigma2_0 = 10000
alpha_0  = 0.01
beta_0   = 0.01

#' Convenience summaries used inside the Gibbs loop.
n    = length(y)
ybar = mean(y)

# --- Initialisation -----------------------------------------------------------

#' Sets the random seed for reproducibility, then randomly initialises
#' mu and tau from broad distributions.  The choice of initial values
#' matters for burn-in length but not for the limiting distribution.
set.seed(13)
mu     = tau = numeric()
sigma2 = 1 / tau
mu[1]     = rnorm(1, 0, 3)
tau[1]    = runif(1, 0, 3)
sigma2[1] = 1 / tau[1]

#' Print the initialised values for inspection.
mu
sqrt(sigma2)

# --- Gibbs sampling loop ------------------------------------------------------

#' Runs Gibbs sampling for the semi-conjugated Normal model.
#' At each iteration:
#'   1. Update mu from its full conditional Normal(mu_1, sigma_1^2)
#'   2. Update tau from its full conditional Gamma(alpha_1, beta_1)
#'   3. Re-express tau on the variance scale
#'
#' The full conditionals are available in closed form because the model is
#' semi-conjugated (Normal-Normal for mu, Gamma-Normal for tau).
nsim = 1000
for (i in 2:nsim) {
  # Full conditional for mu: precision = 1/sigma2_0 + n/sigma2[i-1]
  sigma_1 = sqrt(1 / (1 / sigma2_0 + n / sigma2[i - 1]))
  mu_1    = (mu_0 / sigma2_0 + n * ybar / sigma2[i - 1]) * sigma_1^2
  mu[i]   = rnorm(1, mu_1, sigma_1)

  # Full conditional for tau: shape += n/2, rate += sum((y - mu)^2)/2
  alpha_1   = alpha_0 + n / 2
  beta_1    = beta_0 + sum((y - mu[i])^2) / 2
  tau[i]    = rgamma(1, alpha_1, beta_1)
  sigma2[i] = 1 / tau[i]
}

# --- Gibbs sampling visualisation (fig-mcmc) ----------------------------------
#' Reproduces the step-by-step illustration of Gibbs sampling on a bivariate
#' Normal target.  Uses mvtnorm::dmvnorm() to overlay contours of the (known,
#' for illustration purposes) target distribution.
#' NOTE: This code also re-runs the Gibbs loop from scratch (same seed) so
#' that the visualisation objects (sims, toplot, etc.) are self-contained.

set.seed(13)
y = c(1.2697, 7.7637, 2.2532, 3.4557, 4.1776, 6.4320, -3.6623, 7.7567,
      5.9032, 7.2671, -2.3447, 8.0160, 3.5013, 2.8495, 0.6467, 3.2371,
      5.8573, -3.3749, 4.1507, 4.3092, 11.7327, 2.6174, 9.4942, -2.7639,
      -1.5859, 3.6986, 2.4544, -0.3294, 0.2329, 5.2846)
n = length(y)
ybar    = mean(y)
mu_0    = 0;  sigma2_0 = 10000
alpha_0 = 0.01; beta_0 = 0.01

mu = tau = numeric()
sigma2   = 1 / tau
mu[1]    = rnorm(1, 0, 3)
tau[1]   = runif(1, 0, 3)
sigma2[1] = 1 / tau[1]

nsim = 1000
for (i in 2:nsim) {
  sigma_n = sqrt(1 / (1 / sigma2_0 + n / sigma2[i - 1]))
  mu_n    = (mu_0 / sigma2_0 + n * ybar / sigma2[i - 1]) * sigma_n^2
  mu[i]   = rnorm(1, mu_n, sigma_n)
  alpha_n = alpha_0 + n / 2
  beta_n  = beta_0 + sum((y - mu[i])^2) / 2
  tau[i]  = rgamma(1, alpha_n, beta_n)
  sigma2[i] = 1 / tau[i]
}
sigma = sqrt(sigma2)
sims  = tibble(mu = mu, sigma = sigma)

# Build a bivariate Normal approximation to the posterior for contour overlay
require(mvtnorm)
theta = c(mean(mu), mean(sqrt(sigma2)))
s     = c(var(mu), var(sqrt(sigma2)))
rho   = cor(mu, sqrt(sigma2))
ins   = c(-10, 10)
x1    = seq(theta[1] - abs(ins[1]), theta[1] + abs(ins[1]), length.out = 250)
x2    = seq(theta[2] - abs(ins[2]), theta[2] + abs(ins[2]), length.out = 250)
all   = expand.grid(x1, x2)
Sigma = matrix(c(s[1], s[1] * s[2] * rho, s[1] * s[2] * rho, s[2]), nrow = 2)
toplot = cbind(all, prob = mvtnorm::dmvnorm(all, mean = theta, sigma = Sigma))

# Panel 1: target contour + initialisation point
toplot |> ggplot(aes(x = Var1, y = Var2, z = prob)) + geom_contour() +
  xlim(0, 6) + ylim(0, 6) +
  annotate("text", x = 5, y = 5, label = "$p(\\theta_1, \\theta_2 \\mid y)$") +
  labs(title = "Initialisation ('Iteration 0')") +
  xlab("$\\theta_1$") + ylab("$\\theta_2$") +
  geom_point(data = (sims |> slice(1)), aes(x = mu, y = sigma, z = NULL)) +
  annotate("text", x = sims$mu[1], y = sims$sigma[1],
           label = "$\\theta_1^{(0)},\\theta_2^{(0)}$", vjust = -1, hjust = 0.5)

# Panels 2--3: first and second iteration (full conditionals shown as marginal curves)
# (See source .qmd for the full plotting code for iterations 2 and 3,
#  which follow the same structure with updated conditioning values.)

# Panel 4: trajectory for the first 30 iterations
top30 = sims |> slice_head(n = 30)
plt   = top30 |> ggplot(aes(mu, sigma)) +
  xlim(0, 6) + ylim(0, 6) + labs(title = "After 30 iterations") +
  geom_contour(data = toplot, aes(x = Var1, y = Var2, z = prob), col = "grey90") +
  xlab("$\\theta_1$") + ylab("$\\theta_2$")
for (i in 1:29) {
  plt = plt +
    geom_segment(x = mu[i],     xend = mu[i + 1], y = sigma[i],     yend = sigma[i],
                 arrow = arrow(length = unit(0.20, "cm"), type = "closed"),
                 col = "gray80", inherit.aes = FALSE) +
    geom_segment(x = mu[i + 1], xend = mu[i + 1], y = sigma[i],     yend = sigma[i + 1],
                 arrow = arrow(length = unit(0.20, "cm"), type = "closed"),
                 col = "gray80", inherit.aes = FALSE)
}
plt + geom_text(aes(label = 1:30))

# --- Intuition plots: few vs many iterations (fig-intuition-gibbs) ------------

#' Shows how 15 iterations (too few) vs 1000 (enough) differ in their coverage
#' of the target posterior.  ggExtra::ggMarginal() adds marginal histograms.
nval = 15
dat  = tibble(x = mu, y = sigma) |> top_n(nval)
p    = dat |> ggplot(aes(x, y)) + geom_point() + theme_bw() +
  annotate("text", y = rep(0, nval), x = dat$x, label = "X", col = "red") +
  annotate("text", x = rep(0, nval), y = dat$y, label = "X", col = "red") +
  xlab("$\\theta_1$") + ylab("$\\theta_2$") +
  geom_contour(data = cbind(all, prob = mvtnorm::dmvnorm(all, mean = theta, sigma = Sigma)),
               aes(x = Var1, y = Var2, z = prob))
ggExtra::ggMarginal(p, type = "histogram")

nval = 1000
dat  = tibble(x = mu, y = sigma) |> top_n(nval)
p    = dat |> ggplot(aes(x, y)) + geom_point() + theme_bw() +
  annotate("text", y = rep(0, nval), x = dat$x, label = "X", col = "red") +
  annotate("text", x = rep(0, nval), y = dat$y, label = "X", col = "red") +
  xlab("$\\theta_1$") + ylab("$\\theta_2$") +
  geom_contour(data = cbind(all, prob = mvtnorm::dmvnorm(all, mean = theta, sigma = Sigma)),
               aes(x = Var1, y = Var2, z = prob))
ggExtra::ggMarginal(p, type = "histogram")


# ==============================================================================
# SECTION: MCMC diagnostics
# ==============================================================================

# --- Burn-in and convergence (fig-burnin) -------------------------------------

#' Loads the pre-saved MCMC output (two chains) for the convergence figure.
#' The chains are started from very different initial values to make the
#' burn-in period visible in the traceplot.
load("data/bayesian-computation/convergence.RData")

tibble(iteration = 1:500, chain = res[1:500, 1, 10]) |>
  mutate(label = "Chain1") |>
  bind_rows(tibble(iteration = 1:500, chain = res[1:500, 2, 10], label = "Chain2")) |>
  ggplot(aes(x = iteration, y = chain, col = label)) + geom_line() +
  ylab("Simulated values") + xlab("Iteration") +
  theme(legend.title = element_blank(), legend.position = c(0.75, 0.85),
        legend.background = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = 240, linetype = 2) +
  annotate("text", 90,  -Inf, label = "Burn-in",                  vjust = -2.5) +
  annotate("text", 390, -Inf, label = "Sample after convergence", vjust = -2.5) +
  scale_color_manual(values = c("#FF7F0E", "#1F77B4"))

x = res[, , 10]
colnames(x) = c("Chain 1", "Chain 2")

# --- Potential Scale Reduction (R-hat) ----------------------------------------

#' Computes the Brooks-Gelman-Rubin Potential Scale Reduction (PSR / R-hat).
#' R-hat close to 1 indicates that between-chain and within-chain variances are
#' similar, which is the key signature of convergence.
Rhat = function(x) {
  nsims   = nrow(x)
  nchains = ncol(x)
  sigma2.hat = apply(x, 2, var)
  mu.hat     = apply(x, 2, mean)
  W = mean(sigma2.hat)
  B = nsims * var(mu.hat)
  ((((nsims - 1) / nsims) * W + B / nsims) / W) |> sqrt()
}

# R-hat for the full 500 iterations (includes burn-in -- should be >> 1)
Rhat(x)

# R-hat from iteration 250 onwards (post-convergence -- should be ~1)
Rhat(x[250:500, ])

# --- Effective Sample Size (ESS) ----------------------------------------------

#' Approximates the effective sample size (ESS), accounting for autocorrelation.
#' ESS < nominal sample size because successive MCMC draws are correlated.
#' A rule of thumb is ESS > 4000 for reliable 95% interval estimates.
n_eff = function(x) {
  nsims      = nrow(x)
  nchains    = ncol(x)
  sigma2.hat = apply(x, 2, var)
  mu.hat     = apply(x, 2, mean)
  W = mean(sigma2.hat)
  B = nsims * var(mu.hat)
  nchains * nsims * min((1 / B) * (((nsims - 1) / nsims) * W + B / nsims), 1)
}

n_eff_all = n_eff(x) |> round(2)
n_eff_500 = n_eff(x[250:500, ]) |> round(2)

# Monte Carlo Standard Error (MCSE): sd / sqrt(ESS)
# Smaller MCSE => more reliable posterior mean estimate
mcse_all = format(c(x[, 1], x[, 2]) |> sd() / sqrt(n_eff(x)),
                  digits = 5, nsmall = 2)
mcse_500 = format(c(x[250:500, 1], x[250:500, 2]) |> sd() / sqrt(n_eff(x[250:500, ])),
                  digits = 5, nsmall = 2)

# --- Autocorrelation plots (fig-autocorrelation) ------------------------------

#' Visualises the autocorrelation function of the MCMC chain.
#' High autocorrelation at large lags (as seen before removing burn-in)
#' means the ESS is much lower than the nominal sample size.
c(x) |> bmhe::acfplot() +
  labs(title = "Autocorrelation for all 1000 iterations")

c(x[250:500, ]) |> bmhe::acfplot() +
  labs(title = "Autocorrelation for the last 500 iterations")


# ==============================================================================
# SECTION: Convergence diagnostics for the semi-conjugated model
# ==============================================================================

#' Reruns the semi-conjugated Normal Gibbs sampler with two chains from
#' deliberately different starting points, then checks convergence via R-hat
#' and traceplots.  Because the model is semi-conjugated, convergence is
#' essentially immediate regardless of the starting point.

set.seed(13)
y = c(1.2697, 7.7637, 2.2532, 3.4557, 4.1776, 6.4320, -3.6623, 7.7567,
      5.9032, 7.2671, -2.3447, 8.0160, 3.5013, 2.8495, 0.6467, 3.2371,
      5.8573, -3.3749, 4.1507, 4.3092, 11.7327, 2.6174, 9.4942, -2.7639,
      -1.5859, 3.6986, 2.4544, -0.3294, 0.2329, 5.2846)
n    = length(y)
ybar = mean(y)
mu_0 = 0;  sigma2_0 = 10000
alpha_0 = 0.01; beta_0 = 0.01

#' Wraps the Gibbs loop into a function so it can be called twice
#' with different initial values.
gibbs = function(nsim = 1000) {
  for (i in 2:nsim) {
    sigma_n = sqrt(1 / (1 / sigma2_0 + n / sigma2[i - 1]))
    mu_n    = (mu_0 / sigma2_0 + n * ybar / sigma2[i - 1]) * sigma_n^2
    mu[i]   = rnorm(1, mu_n, sigma_n)
    alpha_n = alpha_0 + n / 2
    beta_n  = beta_0 + sum((y - mu[i])^2) / 2
    tau[i]  = rgamma(1, alpha_n, beta_n)
    sigma2[i] = 1 / tau[i]
  }
  sigma = sqrt(sigma2)
  list(mu = mu, sigma = sigma)
}

sims = list()

# Chain 1: starting values drawn from broad distributions
mu = tau = numeric()
sigma2   = 1 / tau
mu[1]    = rnorm(1, 0, 3)
tau[1]   = runif(1, 0, 3)
sigma2[1] = 1 / tau[1]
sims[[1]] = gibbs(5000)

# Chain 2: starting values drawn from a different region of the parameter space
mu = tau = numeric()
sigma2   = 1 / tau
mu[1]    = rnorm(1, 10, 3)
tau[1]   = runif(1, 5, 10)
sigma2[1] = 1 / tau[1]
sims[[2]] = gibbs(5000)

# R-hat for mu and sigma across the two chains
Rhat(cbind(sims[[1]]$mu,    sims[[2]]$mu))
Rhat(cbind(sims[[1]]$sigma, sims[[2]]$sigma))

# Traceplots for mu and sigma
sims |> bind_rows() |>
  mutate(chain = rep(c("Chain 1", "Chain 2"), each = length(sims[[1]]$mu)),
         x     = rep(1:length(sims[[1]]$mu), 2)) |>
  ggplot(aes(x, mu, col = chain)) + geom_line() + xlab("Iteration") + ylab("") +
  theme(legend.position = c(.3, .1), legend.title = element_blank(),
        legend.background = element_blank()) +
  scale_color_manual(values = c("#FF7F0E", "#1F77B4"))

sims |> bind_rows() |>
  mutate(chain = rep(c("Chain 1", "Chain 2"), each = length(sims[[1]]$mu)),
         x     = rep(1:length(sims[[1]]$mu), 2)) |>
  ggplot(aes(x, sigma, col = chain)) + geom_line() + xlab("Iteration") + ylab("") +
  theme(legend.position = c(.3, .1), legend.title = element_blank(),
        legend.background = element_blank()) +
  scale_color_manual(values = c("#FF7F0E", "#1F77B4"))

# ESS for both parameters
neff_mu    = n_eff(cbind(sims[[1]]$mu,    sims[[2]]$mu))
neff_sigma = n_eff(cbind(sims[[1]]$sigma, sims[[2]]$sigma))


# ==============================================================================
# SECTION: Hamiltonian Monte Carlo (fig-hmc-plot)
# ==============================================================================

#' Illustrative figure only -- no HMC sampler is implemented here.
#' Panel (a): overlays the target density and its negative log-density
#'   (= "potential energy") on a dual y-axis.
#' Panel (b): shows a discretised Hamiltonian trajectory from an initial
#'   position with high potential energy, used to propose a new sample.

dat = tibble(x = seq(-5, 5, .01)) |>
  mutate(d = dnorm(x, 0, 1)) |>
  mutate(nld = -log(d))
scaleFactor = max(dat$d) / max(dat$nld)
steps = tibble(x = seq(-1, 4, .1)) |> mutate(y = -log(dnorm(x)))

# Panel (a): density and potential energy on dual axis
dat |> ggplot(aes(x, d)) + geom_line() + xlab("$\\theta$") + ylab("Density") +
  geom_line(aes(x, nld * scaleFactor), col = "blue") +
  scale_y_continuous(
    "Density",
    sec.axis = sec_axis(~ . / scaleFactor, name = "Potential energy")
  ) +
  theme(axis.title.y.right = element_text(color = "blue")) +
  geom_point(aes(x = 4, y = dnorm(4)), size = 2.5) +
  geom_point(aes(x = 4, y = -log(dnorm(4)) * scaleFactor), size = 2.2) +
  annotate("text", x = 4, y = -log(dnorm(4)) * scaleFactor, size = 3,
           label = "Current position", hjust = 1.1, vjust = -.15)

# Panel (b): Hamiltonian trajectory along the potential energy surface
dat |> ggplot(aes(x, nld)) + ylab("Negative log density") +
  geom_point(data = steps, aes(x = x, y = y), size = 2.2, col = "gray65", shape = 1) +
  geom_point(aes(x = 4, y = -log(dnorm(4))), size = 2.2, col = "gray55") +
  geom_line() + xlab("$\\theta$") +
  geom_point(aes(x = -1, y = -log(dnorm(-1))), size = 2.2) +
  annotate("text", x =  4, y = -log(dnorm(4)),  size = 3, label = "Current position",
           hjust = 1.1, vjust = -.15) +
  annotate("text", x = -1, y = -log(dnorm(-1)), size = 3, label = "Proposal",
           hjust = -.1,  vjust = -.5)
