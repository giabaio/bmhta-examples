#' ---
#' title: "Markov models"
#' desc:  "Bayesian cohort Markov modelling for HTA. Covers: Dirichlet
#'         simplex visualisation; a four-state HIV model combining Multinomial-
#'         Dirichlet transition probability estimation, partial-pooling evidence
#'         synthesis for the relative risk, Weibull survival modelling of
#'         treatment-effect waning, matrix-algebra Markov simulation, discounted
#'         cost and LYG computation; and a three-state cancer Markov model linked
#'         to survival analysis via survHE, with Markov trace and LYG output."
#' ---

library(tidyverse)
library(R2jags)
library(survHE)
library(survHEhmc)
library(ggsimplex)
library(brms)


# ==============================================================================
# SECTION: Dirichlet distribution — simplex visualisations
# ==============================================================================

#' Illustrates how the Dirichlet parameters control mass concentration on
#' the 3-dimensional simplex (equilateral triangle).  Lighter colours indicate
#' higher density; the vertices correspond to the three categories being 100%
#' probable.
#'
#' Uses ggsimplex (for the simplex coordinate system) and brms (for ddirichlet).
#' The custom geom_simplex_canvas() below redraws the simplex canvas with
#' LaTeX-compatible axis labels for a_1, a_2, a_3.

#' Helper: custom simplex canvas geom (vertex labels only; adapted from ggsimplex).
geom_simplex_canvas = function(
    mapping = NULL, data = data.frame(0), stat = "identity",
    position = "identity", na.rm = FALSE, show.legend = NA,
    inherit.aes = FALSE, ...) {
  ggplot2::layer(
    geom        = GeomSimplexCanvas,
    mapping     = mapping,  data = data, stat = stat,
    position    = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params      = list(...)
  )
}
GeomSimplexCanvas = ggplot2::ggproto(
  "GeomSimplexCanvas", ggplot2::GeomPolygon,
  required_aes = c(),
  default_aes  = ggplot2::aes(linewidth = 1, linetype = 1,
                               colour = "black", alpha = 1, fontsize = 18),
  draw_key     = ggplot2::draw_key_polygon,
  draw_panel   = function(data, panel_params, coord) {
    style_params = data
    style_params = coord$transform(style_params, panel_params)
    B      = matrix(c(0, 0, 1, 0, 0.5, sqrt(3) / 2), byrow = TRUE, ncol = 2)
    coords = coord$transform(data.frame(x = B[, 1], y = B[, 2]), panel_params)
    grob   = grid::polygonGrob(
      coords$x, coords$y,
      default.units = "native",
      gp            = grid::gpar(
        col  = scales::alpha(style_params$colour, style_params$alpha),
        fill = scales::alpha("black", 0.0),
        lwd  = style_params$linewidth * .pt,
        lty  = style_params$linetype
      )
    )
    text_style  = grid::gpar(col = "black", fontsize = style_params$fontsize)
    text_grob_1 = grid::textGrob("a1", x = coords$x[1] + 0.021, y = coords$y[1] - 0.018, gp = text_style)
    text_grob_2 = grid::textGrob("a2", x = coords$x[2] - 0.021, y = coords$y[2] - 0.018, gp = text_style)
    text_grob_3 = grid::textGrob("a3", x = coords$x[3],         y = coords$y[3] + 0.026, gp = text_style)
    grid::grobTree(grob, text_grob_1, text_grob_2, text_grob_3)
  }
)

#' Produces one simplex density plot for a given Dirichlet parameter vector.
plot_simplex = function(alpha_vec) {
  df = data.frame(true_model = 1)
  df$Alpha = list(alpha_vec)
  ggplot() +
    coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1)) +
    theme_void() +
    geom_simplex_canvas() +
    ggsimplex::stat_simplex_density(
      data = df, fun = brms::ddirichlet,
      args = alist(Alpha = Alpha)
    )
}

#' Six panels showing different Dirichlet parameterisations:
#'   (1) a=(0.85,0.85,0.85) — diffuse / vague prior
#'   (2) a=(5,5,5)          — symmetric, moderate precision
#'   (3) a=(1,2,3)          — asymmetric, moderate
#'   (4) a=(1,7,2)          — mass pulled towards a_2
#'   (5) a=(10,5,2)         — mass pulled towards a_1
#'   (6) a=(50,50,50)       — very concentrated at barycentre
plot_simplex(c(0.85, 0.85, 0.85))
plot_simplex(c(5,  5,  5))
plot_simplex(c(1,  2,  3))
plot_simplex(c(1,  7,  2))
plot_simplex(c(10, 5,  2))
plot_simplex(c(50, 50, 50))


# ==============================================================================
# SECTION: HIV example — data setup
# ==============================================================================

#' The HIV Markov model (Chancellor et al. 1997) compares:
#'   ZDV      — zidovudine monotherapy
#'   3TC+ZDV  — combination therapy
#'
#' Four health states (ordered by severity):
#'   State A: Compromised CD4  (CD4 in [200; 500])
#'   State B: Low CD4          (CD4 < 200)
#'   State C: AIDS
#'   State D: Death            (absorbing)
#'
#' Only forward transitions are allowed (no recovery).

S      = 4
states = c("Compromised CD4", "Low CD4", "AIDS", "Death")

#' Observed transition counts in the ZDV arm (read by row).
y = matrix(
  c(
    1251, 350, 116, 17,
       0, 731, 512, 15,
       0,   0, 1312, 437,
       0,   0,    0, 469
  ),
  ncol = S, nrow = S, byrow = TRUE
)
colnames(y) = rownames(y) = states

#' Row-wise sample sizes (total transitions observed out of each state)
n = apply(y, 1, sum)
y


# ==============================================================================
# SECTION: Relative risk data (3TC+ZDV vs ZDV)
# ==============================================================================

#' Four published studies comparing 3TC+ZDV to ZDV for AIDS progression.
#' Reported as point estimates and 95% CIs for the RR.
rr  = c(.389, .674, .492, .452)
low = c(.206, .393, .183, .226)
upp = c(.733, 1.156, 1.328, .902)

#' Log-transform the RR to a Normal scale; derive SDs from the 95% CIs.
x  = log(rr)
sd = (log(upp) - log(low)) / (2 * qnorm(.975))


# ==============================================================================
# SECTION: Treatment-effectiveness waning data
# ==============================================================================

#' Pseudo-Binomial data encoding assumptions about how long the 3TC+ZDV
#' treatment remains effective.  r[i] individuals cease to benefit during
#' interval i, out of nr[i] still benefiting at its start.
u     = c(.5, 1, 1.5, 2, 2.5, 3, 4, 5)   # interval endpoints (years)
r     = c(0,  0,   3, 7,  60, 30, 0, 0)   # events (treatment wanes)
I     = length(r)
delta = u - dplyr::lag(u, default = 0)     # interval widths

#' Cumulative at-risk counts: nr[i] = nr[i-1] - r[i-1]
nr = sum(r)
for (i in 2:I) { nr[i] = nr[i - 1] - r[i - 1] }


# ==============================================================================
# SECTION: Cost parameters — Gamma priors for direct care
# ==============================================================================

#' Uses moment-matching to find Gamma(shape, rate) parameters whose 95% interval
#' spans approximately from the "base-case" to the "scenario" cost estimates.
#' bmhe::gammaPar(mean, sd) returns shape and rate.

# Observed cost extremes (base-case and scenario values from Chancellor 1997)
cw  = c(1701, 1774, 6948)   # base-case direct care (£, 1995)
alt = c(2938, 4398, 11223)  # alternative (upper) direct care

pars = list()
pars[[1]] = bmhe::gammaPar(2300, 350)  # Compromised CD4
pars[[2]] = bmhe::gammaPar(2900, 650)  # Low CD4
pars[[3]] = bmhe::gammaPar(9000, 1100) # AIDS

#' Verify the implied distributions cover the required range
bmhe::gammaPar(2300, 350) |>
  (\(x) rgamma(100000, shape = x$shape, rate = x$rate))() |>
  bmhe::stats()

#' Community care: use CV=0.2 (sd = 20% of mean)
bmhe::gammaPar(1055, .2 * 1055)   # Compromised CD4
bmhe::gammaPar(1278, .2 * 1278)   # Low CD4
bmhe::gammaPar(2059, .2 * 2059)   # AIDS

#' Visualise Gamma priors for direct care costs (one panel per state).
for (k in 1:3) {
  xlims = list(c(1000, 4000), c(800, 6000), c(4000, 15000))[[k]]
  print(
    ggplot() +
      stat_function(
        fun  = dgamma,
        args = list(shape = pars[[k]]$shape, rate = pars[[k]]$rate)
      ) +
      xlim(xlims[1], xlims[2]) + xlab("Cost of direct care") +
      geom_segment(
        aes(x = cw[k], xend = alt[k], y = 0, yend = 0),
        col = "blue", linewidth = 1.1
      ) +
      stat_function(
        fun  = dgamma,
        args = list(shape = pars[[k]]$shape, rate = pars[[k]]$rate),
        xlim = c(
          qgamma(.025, shape = pars[[k]]$shape, rate = pars[[k]]$rate),
          qgamma(.975, shape = pars[[k]]$shape, rate = pars[[k]]$rate)
        ),
        geom = "area", fill = "#84CA72", alpha = .2
      ) +
      ylab("")
  )
}


# ==============================================================================
# SECTION: JAGS model — HIV Markov model
# ==============================================================================

#' The model combines four components:
#'
#' 1. Multinomial-Dirichlet model for the ZDV transition matrix Lambda^(1).
#'    Conjugate update: lambda_s | y_s ~ Dirichlet(a_s + y_s).
#'    The Dirichlet prior hyperparameters a are 0 for inadmissible transitions
#'    (e.g. moving from AIDS back to Compromised CD4) and 3 elsewhere.
#'
#' 2. Partial-pooling evidence synthesis for the pooled log RR of progression
#'    under 3TC+ZDV vs ZDV.  Four studies contribute Normally distributed
#'    log-RR estimates x_l with study-specific precisions.  The pooled log RR
#'    is mu.theta; rho = exp(mu.theta) is the pooled RR.
#'
#' 3. Weibull survival model for treatment-effect duration.  Pseudo-Binomial
#'    counts r[i] out of nr[i] encode expert assumptions about when 3TC+ZDV
#'    ceases to be effective.  The Weibull hazard parameterisation:
#'      log(h[i]) = log(alpha) - alpha*log(mu) + (alpha-1)*log(u[i])
#'    The time-varying survival curve surv[j] = exp(-(j/mu)^alpha) is used to
#'    derive the time-varying RR: rho.star[j] = 1 - surv[j]*(1 - rho).
#'    When surv[j] ≈ 1, rho.star ≈ rho (full treatment effect);
#'    when surv[j] ≈ 0, rho.star ≈ 1 (no treatment effect).
#'
#' 4. Gamma priors for the cost of direct care and community care (per state).
#'    Cost of combination therapy: c_drug = c_mono + c_3tc when treatment is
#'    still effective (surv[j] > eps), c_mono otherwise.

model = function() {
  # 1. Multinomial-Dirichlet for ZDV transition matrix
  for (s in 1:S) {
    y[s, 1:S] ~ dmulti(lambda[s, 1:S], n[s])
    lambda[s, 1:S] ~ ddirch(a[s, 1:S])
  }

  # 2. Evidence synthesis for the 3TC+ZDV relative risk
  for (l in 1:L) {
    x[l] ~ dnorm(theta[l], prec.x[l])
    theta[l] ~ dnorm(mu.theta, tau)
  }
  mu.theta    ~ dnorm(0, 0.25)         # pooled log RR (vague prior on logRR scale)
  sigma.theta ~ dexp(2); tau <- pow(sigma.theta, -2)
  rho         <- exp(mu.theta)         # pooled RR

  # 3. Treatment-effectiveness duration (Weibull survival)
  for (i in 1:I) {
    r[i]     ~ dbin(phi[i], nr[i])
    phi[i]   <- 1 - exp(-h[i] * delta[i])       # Pr(wane in interval i)
    log(h[i]) <- log(alpha) - alpha * log(mu) + (alpha - 1) * log(u[i])
  }
  alpha ~ dexp(2)
  mu    ~ dexp(2)

  # Time-varying RR for each Markov cycle j = 0, ..., J
  for (j in 1:(J + 1)) {
    surv[j]      <- exp(-pow(((j - 1) / mu), alpha))
    rho.star[j]  <- 1 - surv[j] * (1 - rho)
  }

  # 4. Costs
  for (s in 1:(S - 1)) {
    c.direct[s] ~ dgamma(shape.dir[s], rate.dir[s])
    c.comm[s]   ~ dgamma(shape.comm[s], rate.comm[s])
  }
  c.direct[S] <- 0; c.comm[S] <- 0   # no cost in Death state

  c.comb[1] <- c.mono + c.3tc
  for (j in 2:(J + 1)) {
    c.comb[j] <- ifelse(surv[(j - 1)] > eps, c.mono + c.3tc, c.mono)
  }
}

#' Data list — all fixed inputs to the JAGS model.
datalist = list(
  y        = y,     S = S,   n = n,
  x        = x,     prec.x = 1 / sd^2,   L = length(x),
  J        = 20,
  a        = matrix(c(3,3,3,3, 0,3,3,3, 0,0,3,3, 0,0,0,3), nrow = S, byrow = TRUE),
  shape.dir  = c(43.184, 19.905, 66.942),  rate.dir  = c(0.019, 0.007, 0.007),
  shape.comm = c(25, 25, 25),              rate.comm = c(0.027, 0.0196, 0.0121),
  u        = c(.5, 1, 1.5, 2, 2.5, 3, 4, 5),
  r        = c(0, 0, 3, 7, 60, 30, 0, 0),
  I        = 8,
  nr       = c(100, 100, 100, 97, 90, 30, 0, 0),
  delta    = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1),
  c.mono   = 2278,   c.3tc = 2086,   eps = 0.001
)

#' Run JAGS.
hiv = jags(
  data               = datalist,
  parameters.to.save = c("lambda", "c.direct", "rho.star", "c.comb", "c.comm", "surv"),
  model.file         = model,
  n.chains           = 2,
  n.iter             = 10000,
  n.burnin           = 5000
)


# ==============================================================================
# SECTION: Post-processing — transition matrix visualisation
# ==============================================================================

#' Applies the time-varying RR to the ZDV transition probabilities to construct
#' the time-varying transition matrix for 3TC+ZDV (Lambda^(2)_j).
#'
#' The rescaling uses the logit-scale formulation (eq-lambda2mat2) for numerical
#' stability:
#'   logit(lambda^(2)) = logit(lambda^(1)) + log(RR) + log(1-lambda^(1)) - log(1-lambda^(1)*RR)
#' then apply the diagonal correction so each row sums to 1.

J            = 20
states.short = c("Compromised", "Low", "AIDS", "Death")
lambda1 = hiv$BUGSoutput$sims.matrix |> as_tibble() |> select(contains("lambda"))
lambda2 = lambda1

#' Build the time-varying Lambda^(2) for visualisation (using dplyr, per cycle).
tmp = lapply(
  1:(J + 1), function(j) {
    lambda2 |>
      mutate(across(
        c(`lambda[1,2]`, `lambda[1,3]`, `lambda[1,4]`,
          `lambda[2,3]`, `lambda[2,4]`, `lambda[3,4]`),
        ~ bmhe::ilogit(
          bmhe::logit(.x) + log(hiv$BUGSoutput$sims.list$rho.star[, j]) +
            log(1 - .x) - log(1 - .x * (1 - hiv$BUGSoutput$sims.list$rho.star[, j]))
        )
      )) |>
      mutate(
        `lambda[1,1]` = 1 - (`lambda[1,2]` + `lambda[1,3]` + `lambda[1,4]`),
        `lambda[2,2]` = 1 - (`lambda[2,3]` + `lambda[2,4]`),
        `lambda[3,3]` = 1 - `lambda[3,4]`
      ) |>
      bmhe::stats() |> as_tibble() |>
      mutate(time = j, From = rep(1:S, S), To = rep(1:S, each = S))
  }
) |>
  bind_rows() |>
  mutate(
    parameter  = paste0("lambda[", From, ",", To, "]"), arm = "3TC+ZDV",
    comp       = paste0(From, "-", To),
    from_lab   = states.short[From], to_lab = states.short[To],
    lab        = as.factor(paste0(from_lab, " -> ", to_lab))
  ) |>
  group_by(comp) |> mutate(id = cur_group_id()) |> ungroup()

#' Plot time-varying transition probabilities for the combination arm.
tmp |>
  dplyr::filter(comp %in% c("1-1","1-2","1-3","1-4","2-2","2-3","2-4","3-3","3-4")) |>
  ggplot(aes(time, mean)) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = "gray80") +
  geom_line() +
  facet_wrap(~fct_reorder(lab, id, .desc = FALSE), ncol = 3, nrow = 3) +
  ylab("Transition probabilities") +
  xlab("Time (Markov cycle)") +
  theme(axis.text.x = element_text(size = 6), strip.text = element_text(size = 7))

#' Posterior summaries: time-varying RR and treatment-effectiveness survival.
bmhe::coefplot(hiv, parameter = "rho.star")$data |>
  mutate(num = row_number(), var = paste0("rho*[", num - 1, "]")) |>
  ggplot(aes(mean, fct_reorder(var, num, .desc = FALSE))) +
  geom_point(size = 1.4) +
  geom_linerange(aes(xmin = low, xmax = upp)) +
  geom_vline(xintercept = 1, col = "gray50") +
  labs(x = "Relative risk") + ylab("Parameter") + xlim(0, 1.1)

bmhe::coefplot(hiv, parameter = "surv")$data |>
  mutate(num = row_number(), j = num - 1) |>
  ggplot(aes(j, mean)) +
  geom_point(size = 1.4) +
  geom_linerange(aes(ymin = low, ymax = upp)) +
  geom_line(
    data = tibble(x = seq(0, 20, .1),
                  y = 1 - pweibull(seq(0, 20, .1), shape = 5.576, scale = 2.707)),
    aes(x, y), col = "gray80"
  ) +
  geom_hline(yintercept = datalist$eps, linetype = 2, col = "grey60") +
  xlab("Time (Markov cycle)") +
  ylab("Probability of treatment effectiveness")


# ==============================================================================
# SECTION: Running the Markov model (matrix algebra)
# ==============================================================================

#' Applies eq-occupancy2 for each MCMC draw: m_{j+1} = m_j * Lambda_j.
#' The arrays m1 and m2 have dimensions [n.sims, J+1, S].
#' R index j=1 corresponds to Markov cycle j=0 (R does not allow 0-indexing).

J     = 20
start = c(1000, 0, 0, 0)   # all 1000 virtual individuals start in State A

m1 = m2 = array(
  0, c(hiv$BUGSoutput$n.sims, (J + 1), S),
  dimnames = list(1:hiv$BUGSoutput$n.sims, paste0("j=", 0:J), states)
)
for (s in 1:S) {
  m1[, 1, s] = start[s]
  m2[, 1, s] = start[s]
}

#' Extract transition matrix simulations and build the 4D lambda2 array.
lambda1 = hiv$BUGSoutput$sims.list$lambda
lambda2 = array(0, dim = c(dim(lambda1), (J + 1)))

for (i in 1:hiv$BUGSoutput$n.sims) {
  for (j in 1:(J + 1)) {
    #' Apply time-varying RR on the logit scale (eq-lambda2mat2)
    lambda2[i, , , j] = bmhe::ilogit(
      bmhe::logit(lambda1[i, , ]) +
        log(hiv$BUGSoutput$sims.list$rho.star[i, j]) +
        log(1 - lambda1[i, , ]) -
        log(1 - lambda1[i, , ] * (1 - hiv$BUGSoutput$sims.list$rho.star[i, j]))
    )
    #' Rescale diagonal so each row sums to 1
    diag(lambda2[i, , , j]) = diag(lambda2[i, , , j]) + (1 - rowSums(lambda2[i, , , j]))
    lambda2[i, S, S, j]     = 1   # Death remains absorbing
  }
}

#' Sanity check: no negative transition probabilities should exist.
which(lambda2 < 0, arr.ind = TRUE)   # expected: empty

#' Run the Markov model: fill m1 and m2 via matrix multiplication.
for (i in 1:hiv$BUGSoutput$n.sims) {
  for (j in 2:(J + 1)) {
    m1[i, j, ] = m1[i, j - 1, ] %*% lambda1[i, , ]
    m2[i, j, ] = m2[i, j - 1, ] %*% lambda2[i, , , j - 1]
  }
}

cbind(m1[1, , ], m2[1, , ]) |> head() |> round(2)

#' Markov trace: mean and 95% interval of state occupancy over time.
m1 |> apply(c(3, 2), bmhe::stats) |> apply(1, t) |> as_tibble() |>
  mutate(state = rep(states, each = J + 1), t = rep(0:J, S),
         intervention = "ZDV", s = rep(1:4, each = J + 1)) |>
  bind_rows(
    m2 |> apply(c(3, 2), bmhe::stats) |> apply(1, t) |> as_tibble() |>
      mutate(state = rep(states, each = J + 1), t = rep(0:J, S),
             intervention = "3TC+ZDV", s = rep(1:4, each = J + 1))
  ) |>
  ggplot(aes(t, mean, col = intervention)) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = intervention), alpha = .25) +
  geom_line(linewidth = .75) + ylab("Individuals") + xlab("Time (Markov cycles)") +
  facet_wrap(~fct_reorder(state, s, .desc = FALSE), nrow = 2, ncol = 2) +
  theme(legend.position = "bottom") + labs(col = "", fill = "") +
  scale_color_manual(values = c("#FF7F0E", "#1F77B4"))


# ==============================================================================
# SECTION: Life years gained (LYG) and incremental effectiveness
# ==============================================================================

#' Occupancy probability arrays (proportions, not counts)
pi1 = m1 / 1000
pi2 = m2 / 1000

#' Length of stay (LOS) in each state per MCMC draw — sum over Markov cycles.
los1 = lapply(1:hiv$BUGSoutput$n.sims, function(i) pi1[i, , ] |> apply(2, sum)) |>
  bind_rows()
los2 = lapply(1:hiv$BUGSoutput$n.sims, function(i) pi2[i, , ] |> apply(2, sum)) |>
  bind_rows()

los1 |> bmhe::stats()
los2 |> bmhe::stats()

#' Incremental effectiveness: LYG = sum of years in all non-Death states.
los1 |> select(-Death) |>
  mutate(LYG1 = rowSums(across(everything()))) |> select(LYG1) |>
  bind_cols(
    los2 |> bind_rows() |> select(-Death) |>
      mutate(LYG2 = rowSums(across(everything()))) |> select(LYG2)
  ) |>
  mutate(Delta.e = LYG2 - LYG1) |>
  bmhe::stats()


# ==============================================================================
# SECTION: Discounted costs
# ==============================================================================

#' Builds the discounted cost tibble for each treatment arm.
#' For each MCMC draw i, we compute the occupancy × (drug + community + direct)
#' cost product for each state and Markov cycle, then apply discount factor.
#'
#' Convention: all cost columns begin with "c_" so mutate(across(starts_with("c_"),
#' ~.x/(1+d)^j)) applies discounting in one line.

d = 0.06   # 6% annual discount rate (Chancellor et al. 1997)

c1 = lapply(1:hiv$BUGSoutput$n.sims, function(i) {
  pi1[i, , ] |> as.data.frame() |>
    cbind(
      j       = 0:J,
      c_drug  = rep(2278, J + 1),   # ZDV cost is constant
      matrix(
        rep(hiv$BUGSoutput$sims.list$c.comm[i, ], J + 1),
        J + 1, S, byrow = TRUE,
        dimnames = list(0:J, paste0("c_comm_", states))
      ),
      matrix(
        rep(hiv$BUGSoutput$sims.list$c.direct[i, ], J + 1),
        J + 1, S, byrow = TRUE,
        dimnames = list(0:J, paste0("c_direct_", states))
      )
    ) |> as_tibble()
}) |> bind_rows() |>
  mutate(across(starts_with("c_"), ~ .x / (1 + d)^j)) |>
  mutate(
    `Compromised CD4` = `Compromised CD4` *
      (`c_drug` + `c_comm_Compromised CD4` + `c_direct_Compromised CD4`),
    `Low CD4`         = `Low CD4` *
      (`c_drug` + `c_comm_Low CD4` + `c_direct_Low CD4`),
    `AIDS`            = `AIDS` *
      (`c_drug` + `c_comm_AIDS` + `c_direct_AIDS`),
    sim               = rep(1:hiv$BUGSoutput$n.sims, each = J + 1)
  ) |>
  select(sim, j, `Compromised CD4`, `Low CD4`, `AIDS`)

#' Combination therapy arm: drug cost varies over the Markov cycles.
c2 = lapply(1:hiv$BUGSoutput$n.sims, function(i) {
  (m2[i, , ] / 1000) |> as.data.frame() |>
    cbind(
      j      = 0:J,
      c_drug = cbind(hiv$BUGSoutput$sims.list$c.comb[i, ]),  # time-varying
      matrix(
        rep(hiv$BUGSoutput$sims.list$c.comm[i, ], J + 1),
        J + 1, S, byrow = TRUE,
        dimnames = list(0:J, paste0("c_comm_", states))
      ),
      matrix(
        rep(hiv$BUGSoutput$sims.list$c.direct[i, ], J + 1),
        J + 1, S, byrow = TRUE,
        dimnames = list(0:J, paste0("c_direct_", states))
      )
    ) |> as_tibble()
}) |> bind_rows() |>
  mutate(across(starts_with("c_"), ~ .x / (1 + d)^j)) |>
  mutate(
    `Compromised CD4` = `Compromised CD4` *
      (`c_drug` + `c_comm_Compromised CD4` + `c_direct_Compromised CD4`),
    `Low CD4`         = `Low CD4` *
      (`c_drug` + `c_comm_Low CD4` + `c_direct_Low CD4`),
    `AIDS`            = `AIDS` *
      (`c_drug` + `c_comm_AIDS` + `c_direct_AIDS`),
    sim               = rep(1:hiv$BUGSoutput$n.sims, each = J + 1)
  ) |>
  select(sim, j, `Compromised CD4`, `Low CD4`, `AIDS`)

#' Cost summaries per state (summed across all Markov cycles).
c1 |> group_by(sim) |>
  summarise(across(`Compromised CD4`:`AIDS`, ~ sum(.x))) |>
  select(-sim) |> bmhe::stats()

c2 |> group_by(sim) |>
  summarise(across(`Compromised CD4`:`AIDS`, ~ sum(.x))) |>
  select(-sim) |> bmhe::stats()

#' Incremental cost Delta_c.
c1 |> mutate(tot.cost = rowSums(across(`Compromised CD4`:`AIDS`))) |>
  select(sim, tot.cost, j) |> group_by(sim) |> summarise(ZDV = sum(tot.cost)) |>
  select(ZDV) |>
  bind_cols(
    c2 |> mutate(tot.cost = rowSums(across(`Compromised CD4`:`AIDS`))) |>
      select(sim, tot.cost, j) |> group_by(sim) |>
      summarise("3TC+ZDV" = sum(tot.cost)) |> select("3TC+ZDV")
  ) |>
  mutate(delta.c = `3TC+ZDV` - ZDV) |>
  bmhe::stats()

#' Cost-over-time plots (annual and cumulative).
c1 |> mutate(tot.cost = rowSums(across(`Compromised CD4`:AIDS))) |>
  select(sim, tot.cost, j) |> group_by(sim) |>
  pivot_wider(names_from = j, names_glue = "j={j}", values_from = tot.cost) |>
  ungroup() |> select(-sim) |> bmhe::stats() |> as_tibble() |>
  mutate(j = 0:J, int = "ZDV") |>
  bind_rows(
    c2 |> mutate(tot.cost = rowSums(across(`Compromised CD4`:AIDS))) |>
      select(sim, tot.cost, j) |> group_by(sim) |>
      pivot_wider(names_from = j, names_glue = "j={j}", values_from = tot.cost) |>
      ungroup() |> select(-sim) |> bmhe::stats() |> as_tibble() |>
      mutate(j = 0:J, int = "3TC+ZDV")
  ) |>
  ggplot(aes(j, mean, col = int)) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = int), alpha = .25) +
  geom_line() + xlab("Time (Markov cycle)") + ylab("Total cost per patient") +
  theme(legend.position = "bottom") + labs(col = "", fill = "") +
  scale_color_manual(values = c("#FF7F0E", "#1F77B4")) +
  scale_fill_manual( values = c("#FF7F0E", "#1F77B4"))


# ==============================================================================
# SECTION: Three-state cancer Markov model (NICE TA174)
# ==============================================================================

#' Loads the TA174 data (digitised NICE appraisal) and restructures it for
#' the three-state cancer Markov model: Pre-progressed -> Progressed -> Death.
#'
#' The three competing risk sets are:
#'   trans==1: Pre-progressed to Progressed  (progression as event)
#'   trans==2: Pre-progressed to Death       (death-without-progression as event)
#'   trans==3: Progressed to Death

data = readRDS("data/markov-models/ta174.rds") |>
  select(patid, treat, prog, death, prog_t, death_t) |>
  mutate(
    treatment = case_when(treat == 0 ~ "FC", TRUE ~ "RFC") |> as.factor()
  ) |>
  rename(id = patid)

#' survHE::make_data_multi_state() restructures the event-history data into
#' the mstate (Andersen-Gill) long format required for multi-transition
#' survival modelling.
risk_set = survHE::make_data_multi_state(data)
risk_set |> as.data.frame() |> head() |> print(digits = 2)

#' Inspect specific risk sets.
risk_set |> dplyr::filter(trans == 1, prog == 0, death == 0) |>
  slice(1:3) |> as.data.frame() |> head() |> print(digits = 2)   # censored in Pre-progression

risk_set |> dplyr::filter(trans == 1, prog == 0, death == 1) |>
  slice(1:3) |> as.data.frame() |> head() |> print(digits = 2)   # die before progression

risk_set |> dplyr::filter(trans == 2, prog == 1, death == 0) |>
  slice(1:3) |> as.data.frame() |> head() |> print(digits = 2)   # progress, not yet dead

#' Fit Bayesian Gompertz survival models for each transition using survHE/rstan.
#' Informative priors on the Gompertz shape are specified to stabilise inference
#' for the sparse Pre-progression-to-Death risk set.
priors = list(gom = list(a_alpha = 1.5, b_alpha = 1.5))

m_12 = fit.models(                         # Pre-progressed -> Progressed
  Surv(time, status) ~ as.factor(treat),
  data   = risk_set |> dplyr::filter(trans == 1),
  distr  = "gom", method = "hmc", priors = priors
)
m_13 = fit.models(                         # Pre-progressed -> Death
  Surv(time, status) ~ as.factor(treat),
  data   = risk_set |> dplyr::filter(trans == 2),
  distr  = "gom", method = "hmc", priors = priors
)
m_23 = fit.models(                         # Progressed -> Death
  Surv(time, status) ~ as.factor(treat),
  data   = risk_set |> dplyr::filter(trans == 3),
  distr  = "gom", method = "hmc", priors = priors
)

#' Survival curve plots for each transition.
plot(m_12, t = seq(0, 120), add.km = TRUE, lab.profile = c("RFC", "FC")) +
  theme(legend.position.inside = c(.6, .75), legend.position = "inside") +
  labs(linetype = "Treatment") + guides(color = "none") + xlim(0, 130)

plot(m_13, t = seq(0, 120), add.km = TRUE, lab.profile = c("RFC", "FC")) +
  theme(legend.position.inside = c(.25, .5), legend.position = "inside") +
  labs(linetype = "Treatment") + guides(color = "none") + xlim(0, 130)

plot(m_23, t = seq(0, 120), add.km = TRUE, lab.profile = c("RFC", "FC")) +
  theme(legend.position.inside = c(.6, .75), legend.position = "inside") +
  labs(linetype = "Treatment") + guides(color = "none") + xlim(0, 130)

#' three_state_mm() runs the three-state Markov model for nsim posterior draws.
#' t defines the time grid (in months); the output mm$m is a tidy tibble.
mm = three_state_mm(m_12, m_13, m_23, t = seq(0, 120), nsim = 1000)

#' Rename treatment profiles for readability.
mm$m = mm$m |>
  mutate(profile = case_when(
    profile == "as.factor(treat)1=0" ~ "FC", TRUE ~ "RFC"
  ))

#' Markov trace: state occupancy over the time horizon.
markov_trace(mm, interventions = c("FC", "RFC")) +
  theme(legend.position = "bottom") +
  scale_fill_manual(
    values = c("#FF7F0E", "#2CA02C", "#1F77B4") |> rev(),
    labels = c("Death", "Progression", "Pre-progression") |> rev(),
    breaks = c(3, 2, 1)
  ) +
  xlab("Markov cycle (months)")

#' Life years gained: sum of months in non-Death states, converted to years.
(mm$m |>
  group_by(profile, sim_idx) |>
  summarise(los = sum(`Pre-progressed`) / 1000 + sum(Progressed) / 1000)) |>
  group_by(profile) |>
  summarise(
    mean   = mean(los) / 12,
    sd     = sd(los) / 12,
    `2.5%` = quantile(los, .025) / 12,
    `97.5%`= quantile(los, .975) / 12
  )
