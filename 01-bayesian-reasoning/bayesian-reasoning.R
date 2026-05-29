#' ---
#' title: "Introduction to Bayesian reasoning"
#' desc:  "Introduces the Bayesian modelling framework through worked examples.
#'         Covers: significance testing and the p-value; Bayes' theorem applied
#'         to Covid-19 testing; prior elicitation for the drug effectiveness
#'         example (Beta and logit-Normal distributions); convergence of expert
#'         opinion with accumulating data; Monte Carlo simulation and the
#'         natural/original parameter scale distinction; the logit-Normal density
#'         and its equivalence with the Beta; and prior predictive (forward)
#'         sampling."
#' ---

library(tidyverse)


# ==============================================================================
# SECTION: Significance testing — the p-value (frequentist approach)
# ==============================================================================

#' Illustrates the deductive/frequentist perspective through a two-sample test.
#' Under the null H0: delta = mu1 - mu0 = 0, the test statistic T = D/sqrt(S^2_D)
#' follows a t distribution.  The p-value is the tail-area probability
#' Pr(T > t_obs | H0), computed using pt().
#'
#' This is contrasted with the Bayesian inductive approach below, which fixes
#' the observed data and reasons about the plausibility of hypotheses.

set.seed(120)
n0 = 80; n1 = 78; n = n0 + n1
mu0 = 25; mu1 = 32
sigma0 = 17; sigma1 = 20

y0 = rnorm(n0, mu0, sigma0)
y1 = rnorm(n1, mu1, sigma1)

y0bar = mean(y0)
y1bar = mean(y1)

#' Sample variances (unbiased estimator, denominator n-1)
s2_0 = sum((y0 - y0bar)^2) / (n0 - 1)   # equivalent to var(y0)
s2_1 = sum((y1 - y1bar)^2) / (n1 - 1)

#' Observed difference and its variance estimate
d    = y1bar - y0bar
s2_D = (s2_0 / n0) + (s2_1 / n1)

#' Test statistic under H0 (delta = 0)
delta = 0
t     = (d - delta) / sqrt(s2_D)

#' One-sided p-value: Pr(T > t_obs | H0)
#' A small p-value indicates the observed result is unlikely under H0.
pt(q = t, df = n - 1, lower.tail = FALSE)

#' Visualise the t sampling distribution and the p-value region.
ggplot(NULL, aes(c(-4, 4))) +
  geom_area(
    stat = "function", fun = dt, args = list(df = n - 1),
    fill = "white", xlim = c(-4, t), col = "black", alpha = 0.15
  ) +
  geom_area(
    stat = "function", fun = dt, args = list(df = n - 1),
    fill = "grey", xlim = c(t, 4), col = "black"
  ) +
  theme_bw() +
  xlab("T") + ylab("Sampling distribution for T") +
  annotate("text", x = t, y = -.005,
           label = paste0("t = ", format(t, digits = 3, nsmall = 2)))


# ==============================================================================
# SECTION: Bayes' theorem — Covid-19 testing
# ==============================================================================

#' Demonstrates how the prior probability of disease (prevalence) modulates the
#' posterior interpretation of a negative test result.
#'
#' Test characteristics (assumed):
#'   Sensitivity:   Pr(+ve | Disease)     = 0.96  =>  Pr(-ve | Disease) = 0.04
#'   Specificity:   Pr(-ve | No disease)  = 0.95
#'
#' Posterior probability of disease given a negative test:
#'   Pr(Disease | -ve) = Pr(Disease) * Pr(-ve | Disease) / Pr(-ve)
#'
#' Even if a p-value analogy would say the negative test rules out disease
#' (p = 0.04), the posterior depends heavily on the prior (prevalence theta).

theta = seq(0, 1, .01)
data  = tibble(
  prior = theta,
  post  = theta * 0.04 / (theta * 0.04 + (1 - theta) * .95)
)

data |>
  ggplot() + geom_line(aes(prior, post)) +
  labs(
    x = "Prior probability of disease, theta",
    y = "Posterior given -ve test"
  ) +
  annotate("text", .1, data |> dplyr::filter(prior == 0.1) |> pull(post),
           label = format(data |> dplyr::filter(prior == 0.1) |> pull(post), digits = 4),
           vjust = -1.2) +
  annotate("text", .4, data |> dplyr::filter(prior == 0.4) |> pull(post),
           label = format(data |> dplyr::filter(prior == 0.4) |> pull(post), digits = 4),
           vjust = -1.2) +
  annotate("text", .8, data |> dplyr::filter(prior == 0.8) |> pull(post),
           label = format(data |> dplyr::filter(prior == 0.8) |> pull(post), digits = 4),
           vjust = -1.2) +
  geom_segment(
    aes(x = .1, y = -Inf, xend = .1,
        yend = data |> dplyr::filter(prior == 0.1) |> pull(post)),
    linetype = "dashed"
  ) +
  geom_segment(
    aes(x = .4, y = -Inf, xend = .4,
        yend = data |> dplyr::filter(prior == 0.4) |> pull(post)),
    linetype = "dashed"
  ) +
  geom_segment(
    aes(x = .8, y = -Inf, xend = .8,
        yend = data |> dplyr::filter(prior == 0.8) |> pull(post)),
    linetype = "dashed"
  )

#' Key insight: if prevalence is 10%, a negative test leaves you with only
#' ~0.46% probability of having the disease — very different from the
#' "p-value" of 0.04 from the frequentist perspective, which ignores the prior.


# ==============================================================================
# SECTION: Two experts — convergence of prior opinions with data
# ==============================================================================

#' Shows how two experts starting from very different priors converge on the
#' same posterior as sample size grows.
#'
#' Expert E1 has a low-effectiveness prior: Beta(3, 12)  => mean ~0.2
#' Expert E2 has a high-effectiveness prior: Beta(12, 3) => mean ~0.8
#'
#' The Beta(a, b) prior is conjugate to the Binomial: after observing y
#' successes in n trials, the posterior is Beta(a+y, b+n-y).
#'
#' With small n (n=10), the two posteriors remain distinct.
#' With large n (n=200), the posteriors effectively coincide.

set.seed(1234)

# Priors
a0 = 3;  b0 = 12   # Expert E1: low effectiveness
c0 = 12; d0 = 3    # Expert E2: high effectiveness

# Small dataset (n=10)
n  = 10
y  = sum(rbinom(n, 1, 0.93))   # simulate y successes
a1 = a0 + y;  b1 = b0 + n - y
c1 = c0 + y;  d1 = d0 + n - y

# Large dataset (n=200)
n1 = 200
y1 = sum(rbinom(n1, 1, 0.93))
a2 = a0 + y1; b2 = b0 + n1 - y1
c2 = c0 + y1; d2 = d0 + n1 - y1

#' Construct a tibble of 95% interval estimates for all prior/posterior stages.
ints = tibble(
  lower = c(
    rbeta(100000, a0, b0) |> quantile(0.025),
    rbeta(100000, c0, d0) |> quantile(0.025),
    rbeta(100000, a1, b1) |> quantile(0.025),
    rbeta(100000, c1, d1) |> quantile(0.025),
    rbeta(100000, a2, b2) |> quantile(0.025),
    rbeta(100000, c2, d2) |> quantile(0.025)
  ),
  upper = c(
    rbeta(100000, a0, b0) |> quantile(0.975),
    rbeta(100000, c0, d0) |> quantile(0.975),
    rbeta(100000, a1, b1) |> quantile(0.975),
    rbeta(100000, c1, d1) |> quantile(0.975),
    rbeta(100000, a2, b2) |> quantile(0.975),
    rbeta(100000, c2, d2) |> quantile(0.975)
  ),
  expert = rep(c("E1", "E2"), 3),
  time   = factor(
    rep(
      c(
        "Prior",
        paste0("Posterior with y1=", y,  ", n1=", n),
        paste0("Posterior with y2=", y1, ", n2=", n1)
      ),
      each = 2
    ),
    levels = c(
      paste0("Posterior with y2=", y1, ", n2=", n1),
      paste0("Posterior with y1=", y,  ", n1=", n),
      "Prior"
    )
  )
)

#' Plot: each row is a stage (prior / small posterior / large posterior).
#' Arrows show how each expert's interval shifts as evidence accumulates.
ints |>
  ggplot() + xlim(0, 1) +
  xlab("Interval estimate for theta") + ylab("") +
  geom_segment(
    aes(x = lower, xend = upper, y = time, yend = time, col = expert),
    linewidth = 1.1
  ) +
  # Arrows showing movement of lower bounds from prior to small-data posterior
  geom_segment(data = tibble(prior = ints$lower[1], posterior = ints$lower[3]),
               aes(x = prior, xend = posterior, y = 3, yend = 2),
               arrow = arrow(length = unit(0.20, "cm"), type = "closed"),
               col = "gray70", size = .2) +
  geom_segment(data = tibble(prior = ints$lower[2], posterior = ints$lower[4]),
               aes(x = prior, xend = posterior, y = 3, yend = 2),
               arrow = arrow(length = unit(0.20, "cm"), type = "closed"),
               col = "gray70", size = .2) +
  geom_segment(data = tibble(prior = ints$upper[1], posterior = ints$upper[3]),
               aes(x = prior, xend = posterior, y = 3, yend = 2),
               arrow = arrow(length = unit(0.20, "cm"), type = "closed"),
               col = "gray70", size = .2) +
  geom_segment(data = tibble(prior = ints$upper[2], posterior = ints$upper[4]),
               aes(x = prior, xend = posterior, y = 3, yend = 2),
               arrow = arrow(length = unit(0.20, "cm"), type = "closed"),
               col = "gray70", size = .2) +
  # Arrows from small-data to large-data posterior
  geom_segment(data = tibble(prior = ints$lower[3], posterior = ints$lower[5]),
               aes(x = prior, xend = posterior, y = 2, yend = 1),
               arrow = arrow(length = unit(0.20, "cm"), type = "closed"),
               col = "gray70", size = .2) +
  geom_segment(data = tibble(prior = ints$lower[4], posterior = ints$lower[6]),
               aes(x = prior, xend = posterior, y = 2, yend = 1),
               arrow = arrow(length = unit(0.20, "cm"), type = "closed"),
               col = "gray70", size = .2) +
  geom_segment(data = tibble(prior = ints$upper[3], posterior = ints$upper[5]),
               aes(x = prior, xend = posterior, y = 2, yend = 1),
               arrow = arrow(length = unit(0.20, "cm"), type = "closed"),
               col = "gray70", size = .2) +
  geom_segment(data = tibble(prior = ints$upper[4], posterior = ints$upper[6]),
               aes(x = prior, xend = posterior, y = 2, yend = 1),
               arrow = arrow(length = unit(0.20, "cm"), type = "closed"),
               col = "gray70", size = .2) +
  theme_classic() +
  theme(
    legend.title    = element_blank(),
    legend.position = c(0.15, 0.15),
    legend.background = element_rect(fill = "transparent"),
    axis.line.y     = element_blank(), axis.ticks.y = element_blank()
  ) +
  scale_color_manual(values = c("#E69F00", "#0072B2"))


# ==============================================================================
# SECTION: Monte Carlo simulation
# ==============================================================================

#' Demonstrates that Monte Carlo quantiles converge to the analytic values as
#' the sample size S increases.  The 2.5% and 97.5% quantiles of Normal(0,1)
#' are qnorm(0.025) ≈ -1.96 and qnorm(0.975) ≈ +1.96.

# Basic demonstration: MC quantiles for Normal(0,1)
S = 100000
y = rnorm(S, 0, 1)
quantile(y, c(0.025, 0.975))

#' Convergence plot: how accuracy improves with increasing S.
S_vals = c(10, 15, 20, 30, 50, 75, 100, 200, 300, 500, 750, 1000,
           2000, 5000, 7500, 10000, 50000, 100000, 500000, 1000000)

y_list = lapply(S_vals, function(s) rnorm(s, 0, 1))
low    = lapply(seq_along(S_vals), function(i) quantile(y_list[[i]], 0.025)) |> unlist()
upp    = lapply(seq_along(S_vals), function(i) quantile(y_list[[i]], 0.975)) |> unlist()

tibble(S = S_vals, quantile = low, type = "0.025 quantile") |>
  bind_rows(tibble(S = S_vals, quantile = upp, type = "0.975 quantile")) |>
  ggplot(aes(S, quantile, col = type)) +
  geom_point() + geom_line() +
  geom_hline(yintercept = qnorm(0.025), linetype = 2) +
  geom_hline(yintercept = qnorm(0.975), linetype = 2) +
  xlab("Number of Monte Carlo samples, S") +
  ylab("Monte Carlo estimates") +
  theme(
    legend.background = element_blank(),
    legend.title      = element_blank(),
    legend.position   = c(.75, .5)
  ) +
  scale_x_continuous(
    trans   = "log",
    labels  = scales::comma,
    breaks  = c(10, 100, 1000, 10000, 100000, 1000000)
  ) +
  scale_color_manual(values = c("#E69F00", "#0072B2"))


# ==============================================================================
# SECTION: Natural- vs original-scale parameters (Gamma example)
# ==============================================================================

#' When we think in terms of mean and sd (natural scale) rather than shape
#' and rate (original scale), we can elicit priors on the natural scale and
#' then use Monte Carlo to derive the implied priors on the original scale.
#'
#' For Y ~ Gamma(nu, gamma):
#'   E[Y] = mu  = nu / gamma
#'   Var[Y] = sigma^2 = nu / gamma^2
#' So: gamma = sqrt(mu / sigma^2)  and  nu = mu * gamma
#'
#' Priors:
#'   mu    ~ log-Normal(5.2, 0.2)    [positive, roughly 100--300 on natural scale]
#'   sigma ~ Exponential(0.35)        [PC-style prior for SD]
#'
#' The implied priors for (nu, gamma) are obtained by the change-of-variables
#' formula -- here we use MC rather than deriving the closed-form expression.

library(tibble)
library(ggplot2)

set.seed(3020)
mu    = rlnorm(10000, 5.2, .2)
sigma = rexp(10000, .35)

# Original-scale parameters derived via MC
gamma = sqrt(mu / sigma^2)
nu    = mu * gamma

# Prior for mu
tibble(mu = seq(0, 600)) |>
  mutate(p1 = dlnorm(seq(0, 600), 5.2, .2)) |>
  ggplot(aes(mu, p1)) + geom_line() + xlab("mu") + ylab("")

# Prior for sigma
tibble(sigma = seq(0, 15, .1)) |>
  mutate(p2 = dexp(seq(0, 15, .1), .35)) |>
  ggplot(aes(sigma, p2)) + geom_line() + xlab("sigma") + ylab("")

# Implied prior for nu (original-scale shape)
tibble(nu) |>
  ggplot(aes(nu)) + geom_density() + xlab("nu") + ylab("") + xlim(0, 30000)

# Implied prior for gamma (original-scale rate)
tibble(gamma) |>
  ggplot(aes(gamma)) + geom_density() + xlab("gamma") + ylab("") + xlim(0, 200)


# ==============================================================================
# SECTION: Drug example — prior elicitation (Beta distribution)
# ==============================================================================

#' Maps the background information "effectiveness roughly between 20% and 60%"
#' onto a Beta(9.2, 13.8) prior, using the moment-matching relationships
#' for the Beta distribution:
#'   mu    = a / (a + b)
#'   sigma = sqrt(a*b / ((a+b)^2 * (a+b+1)))
#'
#' Solving for (a, b) given mu=0.4 and sigma=0.1:
#'   a = mu * (mu*(1-mu)/sigma^2 - 1)    => 9.2
#'   b = (1-mu) * (mu*(1-mu)/sigma^2 - 1) => 13.8

# Prior simulation
nsim  = 10000
alpha = 9.2
beta  = 13.8
theta = rbeta(n = nsim, shape1 = alpha, shape2 = beta)

# Summary statistics for the prior
bmhe::stats(theta)

# Visualise the prior as a histogram
tibble(theta) |>
  ggplot(aes(theta)) +
  geom_histogram(col = "black", fill = "grey") +
  xlab("theta")

#' Compute a specific tail-area probability under the prior.
#' Pr(theta > 0.5): the probability that the drug is effective for the majority.
sum(theta > 0.5) / length(theta)


# ==============================================================================
# SECTION: Prior information vs prior distribution
#          (Beta vs logit-Normal equivalence)
# ==============================================================================

#' The same prior information can be encoded in two equivalent ways:
#'   (a) Beta(9.2, 13.8) directly on theta
#'   (b) Normal(-0.405, 0.4137) on the logit scale, back-transformed to theta
#'
#' The logit-Normal density is derived analytically via the change-of-variable
#' rule; see the chapter for the derivation.

#' Custom logit-Normal density function.
#' x      — value in (0,1) at which to evaluate the density
#' mu     — mean of the Normal on the logit scale
#' sigma  — SD of the Normal on the logit scale
dlogitnorm = function(x, mu, sigma) {
  dnorm(log(x / (1 - x)), mu, sigma) * abs(1 / (x * (1 - x)))
}

#' Verify near-equality of the two distributions at theta = 0.4.
dlogitnorm(0.4, -0.405, 0.413)
dbeta(0.4,     9.2,    13.8)

#' Visualise both densities over [0, 1] — they are effectively identical.
tibble(
  x     = seq(0, 1, .001),
  y     = dlogitnorm(seq(0, 1, .001), -.405, .4137),
  model = "logit-Normal(-0.405, 0.4137)"
) |>
  bind_rows(
    tibble(
      x     = seq(0, 1, .001),
      y     = dbeta(seq(0, 1, .001), 9.2, 13.8),
      model = "Beta(9.2, 13.8)"
    )
  ) |>
  ggplot(aes(x, y, color = model, linetype = model)) +
  geom_line(linewidth = 0.75) +
  theme(
    legend.position.inside = c(.75, .75), legend.position = "inside",
    legend.background      = element_rect(fill = "transparent")
  ) +
  labs(colour = "") + xlab("theta") + ylab("p(theta)") +
  scale_color_manual(values    = c("#E69F00", "#0072B2"), name = "") +
  scale_linetype_manual(values = c("solid",   "dashed"),  name = "")

#' Also overlay the Beta density with a histogram of back-transformed logit-Normal
#' samples — the two should look identical.
dat1 = tibble(x = seq(0, 1, .001), p1 = dbeta(seq(0, 1, .001), 9.2, 13.8))

tibble(p2 = bmhe::ilogit(rnorm(10000, bmhe::logit(.4), 0.413))) |>
  ggplot(aes(p2)) +
  geom_histogram(
    aes(y = after_stat(density),
        color = "logit(theta) ~ Normal(-0.405, 0.413)"),
    fill = "grey"
  ) +
  geom_line(
    data = dat1,
    aes(x, p1, col = "theta ~ Beta(9.2, 13.8)"),
    linewidth = 0.8
  ) +
  xlab("theta") + ylab("") +
  scale_color_manual("", values = c("black", "blue")) +
  theme(
    legend.position   = c(0.75, 0.85),
    legend.background = element_rect(fill = "transparent")
  )


#' Numerical check: the two distributions produce near-identical density values.
#' Re-run theta prior for subsequent use.
nsim  = 10000
alpha = 9.2
beta  = 13.8
theta = rbeta(n = nsim, shape1 = alpha, shape2 = beta)


# ==============================================================================
# SECTION: Forward sampling (prior predictive distribution)
# ==============================================================================

#' Propagates the prior uncertainty on theta forward through the Binomial
#' sampling model to generate a prior predictive distribution for y,
#' the number of successes out of n=20 individuals in a hypothetical trial.
#'
#' Each value of theta (drawn from the prior) generates a different trial
#' outcome; the resulting distribution of y reflects both parameter uncertainty
#' and sampling variability simultaneously.

#' Simulate from the Binomial sampling model
y = rbinom(n = nsim, size = 20, prob = theta)

#' Define a critical threshold indicator (at least 15 successes)
P.crit = (y >= 15)

#' Visualise the prior predictive distribution.
#' Lighter bars indicate values at or above the clinical threshold.
tibble(y = y) |>
  mutate(P.crit = y >= 15) |>
  ggplot(aes(y, fill = P.crit)) +
  geom_bar(stat = "count", col = "black") +
  xlab("y") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#1F77B4", "#FF7F0E"))

#' Numerical summary of the joint prior predictive distribution.
#' Columns: P.crit (indicator), y (predicted successes), theta (probability).
sims = cbind(P.crit, y, theta)
bmhe::stats(sims)

#' Key output:
#'   mean(y)      ≈ 8  — expected successes under prior uncertainty
#'   mean(P.crit) ≈ small — prior probability of meeting the threshold
#'   The distribution of theta is unchanged (no data to update it yet)
