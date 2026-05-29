#' ---
#' title: "Aggregated level data and evidence synthesis"
#' desc:  "Bayesian evidence synthesis using aggregated (summary-level) data.
#'         Covers: prior forward sampling for the logit scale; no-pooling,
#'         complete-pooling and partial-pooling (hierarchical) models for the
#'         magnesium meta-analysis; convergence diagnostics; DIC-based model
#'         comparison and model averaging; shrinkage visualisation; and a
#'         multiparameter evidence synthesis for neuraminidase inhibitors
#'         prophylaxis, including cost-effectiveness analysis via BCEA."
#' ---

library(tidyverse)
library(R2jags)

# Suppress verbose R2jags output during model runs
options(r2j.pb = "none", r2j.quiet = TRUE, r2j.print.program = FALSE)


# ==============================================================================
# SECTION: Data — magnesium meta-analysis (S=16 studies)
# ==============================================================================

#' Each row is one trial. r1/r2 are deaths; n1/n2 are total patients.
#' t=1: placebo arm; t=2: magnesium arm.
trials = tibble(
  trial = seq(1, 16),
  trial_name = c(
    "Morton", "Rasmussen", "Smith", "Abraham", "Felstedt", "Shechter",
    "Ceremuzynski", "Bertschat", "Singh", "Pereira", "Shechter 1",
    "Golf", "Thorgersen", "LIMIT-2", "Shechter 2", "ISIS-4"
  ),
  year = c(
    1984, 1986, 1986, 1987, 1988, 1989, 1989, 1989,
    1990, 1990, 1991, 1991, 1991, 1992, 1995, 1995
  ),
  deaths_pl  = c(2, 23, 7, 1, 8, 9, 3, 1, 11, 7, 12, 13, 8, 118, 17, 2103),
  total_pl   = c(36, 135, 200, 46, 148, 56, 23, 21, 75, 27, 80, 33, 122, 1157, 108, 29039),
  deaths_mag = c(1, 9, 2, 1, 10, 1, 1, 0, 6, 1, 2, 5, 4, 90, 4, 2216),
  total_mag  = c(40, 135, 200, 48, 150, 59, 25, 22, 76, 27, 89, 23, 130, 1159, 107, 29011)
)
colnames(trials) = c(
  "ID", "Name", "Year",
  "Deaths (r1)", "Total (n1)", "Deaths (r2)", "Total (n2)"
)

#' Extract vectors for JAGS.
data_mag = as.data.frame(trials)
r1 = data_mag[, 4]; n1 = data_mag[, 5]
r2 = data_mag[, 6]; n2 = data_mag[, 7]

#' Shared data list for all three magnesium models.
data = list(r1 = r1, r2 = r2, n1 = n1, n2 = n2, S = nrow(trials))

#' Initialise the model list.
model = list()


# ==============================================================================
# SECTION: Prior forward sampling — logit scale sanity check
# ==============================================================================

#' Before running any model, check that the chosen priors cover the full [0,1]
#' range for the implied probabilities.  This is a necessary condition for the
#' priors to be "vague" on the probability scale, even if they look tight on
#' the logit scale.

#' Prior for baseline log-odds alpha_s ~ Normal(0, sd=4).
#' Implied prior on the placebo arm probability of death pi_s = logit^{-1}(alpha_s).
tibble(alpha = rnorm(10000, 0, 4)) |>
  mutate(pi = exp(alpha) / (1 + exp(alpha))) |>
  select(pi) |>
  bmhe::stats()

#' Adding delta_s ~ Normal(0, sd=2) gives the implied prior for the magnesium
#' arm: pi_s = logit^{-1}(alpha_s + delta_s).
tibble(alpha = rnorm(10000, 0, 4), delta = rnorm(10000, 0, 2)) |>
  mutate(pi = exp(alpha + delta) / (1 + exp(alpha + delta))) |>
  select(pi) |>
  bmhe::stats()

#' Both distributions effectively span [0,1], confirming the priors are
#' sufficiently vague when mapped back to the probability scale.

#' Visualise the implied priors.
tibble(alpha = rnorm(10000, 0, 4)) |>
  mutate(p = exp(alpha) / (1 + exp(alpha))) |>
  ggplot(aes(p)) +
  geom_histogram(fill = "gray", col = "black") +
  xlab("pi_s  [placebo arm]") + ylab("")

tibble(alpha = rnorm(10000, 0, 4), delta = rnorm(10000, 0, 2)) |>
  mutate(p = exp(alpha + delta) / (1 + exp(alpha + delta))) |>
  ggplot(aes(p)) +
  geom_histogram(fill = "gray", col = "black") +
  xlab("pi_s  [magnesium arm]") + ylab("")

#' Forward sampling check for the Gamma(0.001, 0.001) precision prior,
#' illustrating why it is problematic for hierarchical SDs: almost all the mass
#' sits extremely close to 0 (i.e. the SD is pulled towards very large values).
tau = rgamma(10000, shape = 0.001, rate = 0.001)
sum(tau > 0.001) / length(tau)   # effectively 0 — highly concentrated at tau~0


# ==============================================================================
# SECTION: Model 1 — No-pooling
# ==============================================================================

#' Each study is modelled independently.  alpha_s and delta_s are study-specific
#' and receive no common structure.  This is the most flexible model but cannot
#' produce a single pooled estimate of the magnesium treatment effect.
#'
#' Priors (on the logit scale, verified above):
#'   alpha_s ~ Normal(0, sd=4)   [precision = 0.0625]
#'   delta_s ~ Normal(0, sd=2)   [precision = 0.25]
#'
#' Thinning (n.thin=4) is applied because, with no information sharing, small
#' studies (especially study 8: Bertschat, 1 death total) can mix slowly.

no_pooling = function() {
  for (s in 1:S) {
    r1[s] ~ dbin(pi1[s], n1[s])
    r2[s] ~ dbin(pi2[s], n2[s])
    logit(pi1[s]) <- alpha[s]
    logit(pi2[s]) <- alpha[s] + delta[s]
    alpha[s] ~ dnorm(0, 0.0625)   # sd=4 on logit scale
    delta[s] ~ dnorm(0, 0.25)     # sd=2 on logit scale
  }
}

model$no_pooling = jags(
  data               = data,
  parameters.to.save = c("alpha", "delta"),
  inits = function() {
    list(alpha = rnorm(16, 0, .5), delta = rnorm(16, 0, .5))
  },
  n.chains = 2, n.iter = 10000, n.burnin = 2000, n.thin = 4,
  DIC        = TRUE,
  model.file = no_pooling
)

# --- Diagnostics for no-pooling model ---

#' R-hat plot: flag nodes with R-hat > 1.1.
bmhe::diagplot(model$no_pooling)

#' ESS plot: highlight nodes below 10% of nominal sample size (400).
#' These may need longer chains or different priors.
bmhe::diagplot(model$no_pooling, what = "n.eff", label = TRUE) +
  geom_hline(yintercept = 400, linetype = "dashed") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 400,
           alpha = 0.2, fill = "gray70")

#' Traceplots for nodes with n.eff < 400.
bmhe::traceplot(
  model$no_pooling,
  parameter = which(model$no_pooling$BUGSoutput$summary[, "n.eff"] < 400) |> names()
)


# ==============================================================================
# SECTION: Model 2 — Complete pooling
# ==============================================================================

#' A single common treatment effect d is assumed across all studies.
#' This is the "fixed effects" model; it can produce a pooled OR but is
#' sensitive to heterogeneity and highly influenced by large studies (ISIS-4).
#'
#' OR = exp(d) is monitored directly for easy interpretation.

complete_pooling = function() {
  for (s in 1:S) {
    r1[s] ~ dbin(pi1[s], n1[s])
    r2[s] ~ dbin(pi2[s], n2[s])
    logit(pi1[s]) <- alpha[s]
    logit(pi2[s]) <- alpha[s] + d
    alpha[s] ~ dnorm(0, 0.0625)   # study-specific baseline, not pooled
  }
  d  ~ dnorm(0, 0.25)
  OR <- exp(d)
}

model$complete_pooling = jags(
  data               = data,
  parameters.to.save = c("alpha", "d", "OR"),
  inits = function() {
    list(alpha = rnorm(16, 0, .5), d = rnorm(1))
  },
  n.chains = 2, n.iter = 10000, n.burnin = 2000, n.thin = 4,
  DIC        = TRUE,
  model.file = complete_pooling
)

# --- Diagnostics for complete pooling model ---
bmhe::diagplot(model$complete_pooling)
bmhe::diagplot(model$complete_pooling, what = "n.eff", label = TRUE) +
  geom_hline(yintercept = 400, linetype = "dashed") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 400,
           alpha = 0.2, fill = "gray70")


# ==============================================================================
# SECTION: Model 3 — Partial pooling (hierarchical / random effects)
# ==============================================================================

#' The study-specific log ORs delta_s are drawn from a common Normal(d, sigma_delta).
#' This is the "random effects" model: it borrows strength across studies while
#' still allowing for heterogeneity.
#'
#' Key additions vs complete pooling:
#'   sigma.delta ~ Exponential(2.31)   PC prior: Pr(sigma_delta > 1) ≈ 0.1
#'   delta.pred  ~ Normal(d, prec)     predictive distribution for a new study
#'
#' delta.pred is the best guess for the effect in an *unobserved* future study
#' assumed exchangeable with the 16 analysed; it is centred on d but wider.

partial_pooling = function() {
  for (s in 1:S) {
    r1[s] ~ dbin(pi1[s], n1[s])
    r2[s] ~ dbin(pi2[s], n2[s])
    logit(pi1[s]) <- alpha[s]
    logit(pi2[s]) <- alpha[s] + delta[s]
    alpha[s] ~ dnorm(0, 0.0625)
    delta[s] ~ dnorm(d, prec)    # structured: drawn from common distribution
  }
  d           ~ dnorm(0, 0.25)
  sigma.delta ~ dexp(2.31)       # PC prior: Pr(sigma_delta > 1) ≈ 0.1
  prec        <- pow(sigma.delta, -2)
  OR          <- exp(d)
  delta.pred  ~ dnorm(d, prec)   # predictive log OR for a new study
}

model$partial_pooling = jags(
  data               = data,
  parameters.to.save = c("delta", "d", "delta.pred", "OR", "sigma.delta", "alpha"),
  inits = function() {
    list(
      alpha = rnorm(16), delta = rnorm(16),
      d = rnorm(1), sigma.delta = runif(1)
    )
  },
  n.chains = 2, n.iter = 10000, n.burnin = 2000, n.thin = 4,
  DIC        = TRUE,
  model.file = partial_pooling
)

# --- Diagnostics for partial pooling model ---

#' With information sharing, convergence is much better than for no-pooling.
bmhe::diagplot(model$partial_pooling) +
  geom_text(
    aes(
      label = ifelse(
        model$partial_pooling$BUGSoutput$summary[, "Rhat"] > 1.1,
        model$partial_pooling$BUGSoutput$summary |> rownames(), ""
      )
    ), vjust = -1
  )
bmhe::diagplot(model$partial_pooling, what = "n.eff", label = TRUE) +
  geom_hline(yintercept = 400, linetype = "dashed")


# ==============================================================================
# SECTION: DIC-based model comparison
# ==============================================================================

#' Extracts pD from each model via rjags::dic.samples().
#' pD is the BUGS-style penalty (estimated number of effective parameters),
#' which accounts for correlation induced by the hierarchical structure.
#' For the partial pooling model, pD < nominal parameter count (35) because
#' the delta_s parameters are correlated — this is the "borrowing of information".

pD_no_pooling = rjags::dic.samples(
  model$no_pooling$model, 4000, progress.bar = "none", quiet = TRUE
)
pD_complete_pooling = rjags::dic.samples(
  model$complete_pooling$model, 4000, progress.bar = "none", quiet = TRUE
)
pD_partial_pooling = rjags::dic.samples(
  model$partial_pooling$model, 4000, progress.bar = "none", quiet = TRUE
)

#' Summary table: pV (JAGS default) and DIC for the three models.
tibble(
  Model = c("No-pooling", "Complete pooling", "Partial pooling"),
  pV    = c(
    model$no_pooling$BUGSoutput$pV      |> round(2),
    model$complete_pooling$BUGSoutput$pV |> round(2),
    model$partial_pooling$BUGSoutput$pV  |> round(2)
  ),
  DIC = c(
    model$no_pooling$BUGSoutput$DIC      |> round(2),
    model$complete_pooling$BUGSoutput$DIC |> round(2),
    model$partial_pooling$BUGSoutput$DIC  |> round(2)
  ),
  Dbar = c(
    model$no_pooling$BUGSoutput$summary["deviance", "mean"]      |> round(2),
    model$complete_pooling$BUGSoutput$summary["deviance", "mean"] |> round(2),
    model$partial_pooling$BUGSoutput$summary["deviance", "mean"]  |> round(2)
  ),
  sd_D = c(
    model$no_pooling$BUGSoutput$summary["deviance", "sd"]      |> round(2),
    model$complete_pooling$BUGSoutput$summary["deviance", "sd"] |> round(2),
    model$partial_pooling$BUGSoutput$summary["deviance", "sd"]  |> round(2)
  )
) |> tinytable::tt() |> tinytable::style_tt(j = 2:5, align = "c")

#' Stacked bar chart partitioning DIC into Dbar and penalty, for both pV and pD.
tibble(
  x     = c(
    model$no_pooling$BUGSoutput$summary["deviance", "mean"],
    model$complete_pooling$BUGSoutput$summary["deviance", "mean"],
    model$partial_pooling$BUGSoutput$summary["deviance", "mean"],
    model$no_pooling$BUGSoutput$pV,
    model$complete_pooling$BUGSoutput$pV,
    model$partial_pooling$BUGSoutput$pV
  ),
  Model = rep(c("No-pooling", "Complete pooling", "Partial pooling"), 2),
  type  = rep(c("Dbar", "pV"), each = 3),
  facet = "Using pV as penalty"
) |>
  bind_rows(tibble(
    x = c(
      pD_no_pooling$deviance      |> sum(),
      pD_complete_pooling$deviance |> sum(),
      pD_partial_pooling$deviance  |> sum(),
      pD_no_pooling$penalty       |> sum(),
      pD_complete_pooling$penalty  |> sum(),
      pD_partial_pooling$penalty   |> sum()
    ),
    Model = rep(c("No-pooling", "Complete pooling", "Partial pooling"), 2),
    type  = rep(c("Dbar", "pD"), each = 3),
    facet = "Using pD as penalty"
  )) |>
  ggplot(aes(x = Model, y = x, fill = type)) +
  geom_bar(position = "stack", stat = "identity") +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  ylab("DIC") +
  geom_text(
    aes(label = x |> round(2)), size = 3.5,
    position = position_stack(vjust = .5)
  ) +
  facet_wrap(~facet, scales = "fixed") +
  scale_fill_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C"))


# ==============================================================================
# SECTION: DIC-based model averaging
# ==============================================================================

#' Rather than selecting the single best model, weights proportional to
#' exp(-0.5 * ΔDIC_h) combine all three models.
#' In this case the partial pooling model dominates; the complete pooling
#' model receives a modest non-zero weight because ΔDIC is around 5–7.

dic        = unlist(lapply(model, function(x) x$BUGSoutput$DIC))
delta.dic  = abs(min(dic) - dic)
weights    = exp(-.5 * delta.dic) / sum(exp(-.5 * delta.dic))
format(weights, scientific = FALSE)


# ==============================================================================
# SECTION: Results summary table (complete vs partial pooling)
# ==============================================================================

#' Summarises OR estimates, between-study SD and model fit statistics.
tab = tibble(
  Model = c("Complete pooling", "Partial pooling"),
  estimate = c(
    paste0(
      model$complete_pooling$BUGSoutput$summary["OR", "mean"] |> round(2), " (",
      model$complete_pooling$BUGSoutput$summary["OR", "2.5%"] |> round(2), "; ",
      model$complete_pooling$BUGSoutput$summary["OR", "97.5%"] |> round(2), ")"
    ),
    paste0(
      model$partial_pooling$BUGSoutput$summary["OR", "mean"] |> round(2), " (",
      model$partial_pooling$BUGSoutput$summary["OR", "2.5%"] |> round(2), "; ",
      model$partial_pooling$BUGSoutput$summary["OR", "97.5%"] |> round(2), ")"
    )
  ),
  sigma_delta = c(
    "---",
    paste0(
      model$partial_pooling$BUGSoutput$summary["sigma.delta", "mean"] |> round(2), " (",
      model$partial_pooling$BUGSoutput$summary["sigma.delta", "2.5%"] |> round(2), "; ",
      model$partial_pooling$BUGSoutput$summary["sigma.delta", "97.5%"] |> round(2), ")"
    )
  ),
  deviance = c(
    model$complete_pooling$BUGSoutput$summary["deviance", "mean"] |> round(2),
    model$partial_pooling$BUGSoutput$summary["deviance", "mean"]  |> round(2)
  ),
  pD  = c(
    pD_complete_pooling$penalty |> sum() |> round(2),
    pD_partial_pooling$penalty  |> sum() |> round(2)
  ),
  dic = c(
    (model$complete_pooling$BUGSoutput$summary["deviance", "mean"] +
       pD_complete_pooling$penalty |> sum()) |> round(2),
    (model$partial_pooling$BUGSoutput$summary["deviance", "mean"] +
       pD_partial_pooling$penalty |> sum()) |> round(2)
  )
)
colnames(tab) = c(
  "Model", "OR (95% interval)", "Between-study sd (sigma_delta)",
  "Dbar", "pD", "DIC"
)
tab |> tinytable::tt() |> tinytable::theme_tt("resize", width = 1)


# ==============================================================================
# SECTION: Forest plot — no-pooling vs partial pooling vs complete pooling
# ==============================================================================

#' Overlays posterior means and 95% intervals for the study-specific log ORs
#' (delta_s) from all three models, plus the partial pooling grand mean (d)
#' and predictive distribution (delta.pred).
#'
#' Key visual insight: for small studies, the partial pooling estimate is pulled
#' towards d (shrinkage), while for large studies (ISIS-4) the two models agree.

toplot = as_tibble(model$partial_pooling$BUGSoutput$sims.matrix) |>
  select(contains("delta["), d, delta.pred)

tp = tibble(
  mean = apply(toplot, 2, mean),
  q1   = apply(toplot, 2, quantile, 0.025),
  q2   = apply(toplot, 2, quantile, 0.975)
) |>
  bind_cols(name = c(paste(trials$Name, trials$Year), "posterior mean", "predictive mean")) |>
  select(name, everything()) |>
  mutate(model = "Partial pooling", y = row_number() + .25)

toplot = as_tibble(model$no_pooling$BUGSoutput$sims.matrix) |>
  select(contains("delta["))
tp2 = tibble(
  mean = apply(toplot, 2, mean),
  q1   = apply(toplot, 2, quantile, 0.025),
  q2   = apply(toplot, 2, quantile, 0.975)
) |>
  bind_cols(name = paste(trials$Name, trials$Year)) |>
  select(name, everything()) |>
  mutate(model = "No pooling", y = row_number())

toplot = as_tibble(model$complete_pooling$BUGSoutput$sims.matrix) |> select(d)
tp3 = tibble(
  mean = apply(toplot, 2, mean),
  q1   = apply(toplot, 2, quantile, 0.025),
  q2   = apply(toplot, 2, quantile, 0.975)
) |>
  bind_cols(name = "") |>
  select(name, everything()) |>
  mutate(model = "Complete pooling", y = 17.25)

tp = tp |>
  bind_rows(tp2) |>
  mutate(y = case_when(y > 17 ~ y - .25, TRUE ~ y)) |>
  bind_rows(tp3)

labs = c(paste0("delta[", 1:16, "]"), "d", "delta_pred")

tp |> ggplot() +
  geom_vline(xintercept = 0, size = 0.9, linetype = "dashed") +
  geom_segment(aes(x = q1, xend = q2, y = 19 - y, yend = 19 - y), linetype = "dashed") +
  geom_point(aes(x = mean, y = 19 - y, color = model, shape = model), size = 1.2) +
  theme_bw() +
  labs(x = "log-odds ratio", y = NULL) +
  scale_y_continuous(
    breaks  = 18:1,
    labels  = c(paste(trials$Name, trials$Year), "Posterior mean", "Predictive mean")
  ) +
  theme(text = element_text(size = 10), legend.title = element_blank(),
        legend.position = "bottom") +
  coord_cartesian(xlim = c(-4, 4)) +
  annotate("text", x = 4, y = 18:1, label = labs, size = 3) +
  scale_color_manual(values = c("#FF7F0E", "#2CA02C", "#1F77B4"))


# ==============================================================================
# SECTION: Shrinkage visualisation
# ==============================================================================

#' Plots alpha_s (baseline) vs delta_s (log OR), with arrows from the
#' no-pooling to the partial pooling estimate showing the direction and
#' magnitude of shrinkage towards the grand mean d.
#' Small studies move further; large studies barely move.

tp_ranef = tibble(
  x     = model$partial_pooling$BUGSoutput$summary[
    grep("alpha", rownames(model$partial_pooling$BUGSoutput$summary)), "mean"
  ],
  y     = model$partial_pooling$BUGSoutput$summary[
    grep("delta", rownames(model$partial_pooling$BUGSoutput$summary))[1:16], "mean"
  ],
  model = "Partial pooling"
)
tp_indep = tibble(
  x     = model$no_pooling$BUGSoutput$summary[
    grep("alpha", rownames(model$no_pooling$BUGSoutput$summary)), "mean"
  ],
  y     = model$no_pooling$BUGSoutput$summary[
    grep("delta", rownames(model$no_pooling$BUGSoutput$summary)), "mean"
  ],
  model = "No pooling"
)

tp_all = tp_ranef |> bind_rows(tp_indep)
tp3 = tp_indep |>
  rename(x_indep = x, y_indep = y) |> select(-model) |>
  bind_cols(tp_ranef |> rename(x_ranef = x, y_ranef = y) |> select(-model)) |>
  mutate(
    trial = row_number(),
    diff  = case_when(y_indep > y_ranef ~ 1, y_indep <= y_ranef ~ -1),
    y_ranef2 = y_ranef + diff * 0.007
  )

ggplot() +
  geom_segment(
    data = tp3,
    aes(x = x_indep, y = y_indep, xend = x_ranef, yend = y_ranef2),
    arrow = arrow(length = unit(0.12, "cm"), type = "closed")
  ) +
  geom_point(data = tp_all, aes(x, y, color = model, shape = model), size = 1.2) +
  theme_bw() +
  labs(x = "alpha_s", y = "delta_s") +
  annotate("text", x = tp3$x_indep, y = tp3$y_indep,
           label = tp3$trial, size = 3, hjust = 1.7, vjust = 0) +
  geom_point(
    aes(
      x = mean(tp3$x_indep),
      y = model$partial_pooling$BUGSoutput$summary["d", "mean"]
    ),
    size = 1.2, col = "#FF7F0E"
  ) +
  annotate(
    "text", x = mean(tp3$x_indep),
    y   = model$partial_pooling$BUGSoutput$summary["d", "mean"],
    label = "d", size = 3, vjust = -1.25
  ) +
  scale_color_manual(label = c("No pooling", "Partial pooling"),
                     values = c("#2CA02C", "#1F77B4")) +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  scale_shape_manual(values = c(17, 15))


# ==============================================================================
# SECTION: Influenza — multiparameter evidence synthesis and economic evaluation
# ==============================================================================

#' Combines two evidence streams to estimate the cost-effectiveness of
#' neuraminidase inhibitor (NI) prophylaxis vs status quo:
#'   (1) S=6 head-to-head RCTs of NIs vs placebo -> pooled log OR (mu.delta)
#'   (2) H=9 single-arm studies of influenza incidence -> pooled background
#'       probability p1 = logit^{-1}(mu.beta)
#'
#' The two modules are connected using the log-odds relationship:
#'   logit(p2) = logit(p1) + mu.delta
#'
#' Cost and benefit parameters (c.inf, l) are given informative priors derived
#' from published resource-use data.

# --- Head-to-head study data ---
tab = tibble(
  author = c(
    "Monto et al. (1999)", "Hayden et al. (2000)", "Kaiser et al. (2000)",
    "Hayden et al. (1999)", "Hayden et al. (1999)", "Welliver et al. (2001)"
  ),
  age = c("18-64", "5+", "13-65", "18-65", "18-65", "12-85"),
  n1  = c(554, 423, 144, 268, 251, 462),
  r1  = c(34, 40, 9, 19, 6, 34),
  n2  = c(553, 414, 144, 268, 252, 493),
  r2  = c(11, 7, 3, 3, 3, 4)
)

# --- Single-arm incidence data ---
y = c(0, 6, 5, 6, 25, 18, 14, 3, 27)   # influenza cases
m = c(23, 241, 159, 137, 519, 298, 137, 24, 132)  # totals
H = length(y)

# --- Fixed cost parameters ---
unit.cost.drug = 2.40     # daily unit cost of NIs (£)
length.treat   = 6 * 7   # 6-week course (days)
c.gp           = 19       # cost of GP visit (£)
vat            = 1.175    # VAT multiplier
c.ni           = unit.cost.drug * length.treat * vat  # total NI cost

# --- Informative prior parameters for the influenza treatment cost ---
mu.inf    = 16.78          # mean cost of treating influenza (£)
sigma.inf = 2.34           # SD
tau.inf   = 1 / sigma.inf^2  # JAGS precision (= 1 / 2.34^2 ≈ 0.1826)

# --- Log-Normal parameters for duration of infection ---
#' bmhe::lognPar() maps (mean, sd) on the natural scale to (mulog, sigmalog).
#' Target: mean = 8.2 days, sd = 2 days (sqrt(2) in original paper).
mu.l    = bmhe::lognPar(8.2, 2)$mulog
sigma.l = bmhe::lognPar(8.2, 2)$sigmalog
#' Results: mu.l ≈ 2.075, sigma.l ≈ 0.2403 -> precision ≈ 17.3051

#' Combined data list for all three influenza models.
dataflu = list(
  r1 = tab$r1, r2 = tab$r2, n1 = tab$n1, n2 = tab$n2, S = nrow(tab),
  y  = y, m = m, H = H
)

# --- Base-case model: alpha_s and beta_h modelled separately ---

#' alpha_s (head-to-head study baselines) are given independent Normal(0, sd=10)
#' priors — effectively no-pooling for the baselines.
#' beta_h (single-arm incidence studies) are pooled via a shared Normal(mu.beta, tau.beta).
#' delta_s are partially pooled via Normal(mu.delta, tau.delta).
#'
#' Probabilities:
#'   p1 = logit^{-1}(mu.beta)           [baseline probability of influenza]
#'   p2 = logit^{-1}(logit(p1) + mu.delta)  [probability under NIs]
#'
#' Cost and benefit parameters are modelled with informative priors:
#'   c.inf ~ Normal(16.78, precision=0.1826)   [NICE unit cost of influenza]
#'   l     ~ log-Normal(2.075, precision=17.3051)  [duration of infection, days]

influenza = function() {
  # 1. Relative effectiveness of NIs (S head-to-head studies)
  for (s in 1:S) {
    r1[s] ~ dbin(pi1[s], n1[s])
    logit(pi1[s]) <- alpha[s]
    r2[s] ~ dbin(pi2[s], n2[s])
    logit(pi2[s]) <- alpha[s] + delta[s]
    alpha[s] ~ dnorm(0, 0.1)              # independent baseline (no-pooling)
    delta[s] ~ dnorm(mu.delta, tau.delta) # partially pooled log OR
  }
  mu.delta    ~ dnorm(0, 0.1)
  sigma.delta ~ dexp(1.65)                # PC prior: Pr(sigma_delta > 1) ≈ 0.2
  tau.delta   <- pow(sigma.delta, -2)
  or          <- exp(mu.delta)

  # 2. Incidence of influenza in the general population (H single-arm studies)
  for (h in 1:H) {
    y[h] ~ dbin(theta[h], m[h])
    logit(theta[h]) <- beta[h]
    beta[h] ~ dnorm(mu.beta, tau.beta)    # partially pooled baseline risk
  }
  mu.beta    ~ dnorm(0, 0.1)
  sigma.beta ~ dexp(1.65)
  tau.beta   <- pow(sigma.beta, -2)

  # 3. Combine to produce arm-specific probabilities
  logit(p1) <- mu.beta                    # baseline probability without NIs
  logit(p2) <- logit(p1) + mu.delta       # probability with NIs

  # 4. Cost and duration parameters (informative priors)
  c.inf ~ dnorm(16.78, 0.1826)            # cost per influenza episode (£)
  l     ~ dlnorm(2.075, 17.3051)          # duration of infection (days)
}

flu = jags(
  data               = dataflu,
  parameters.to.save = c("p1", "p2", "or", "c.inf", "l"),
  inits = function() {
    list(alpha = rnorm(dataflu$S, 0, .5), delta = rnorm(dataflu$S, 0, .5))
  },
  n.chains = 2, n.iter = 25000, n.burnin = 5000, n.thin = 10,
  DIC        = TRUE, pD = TRUE,
  model.file = influenza
)
print(flu)

# --- Variant 1: pool all baseline risks (alpha_s and beta_h exchangeable) ---

#' Relaxes the assumption that head-to-head and single-arm baselines are
#' independent by modelling all of them from Normal(mu.beta, tau.beta).
#' This is appropriate if the two sets of studies are sufficiently similar.

influenza2 = function() {
  for (s in 1:S) {
    r1[s] ~ dbin(pi1[s], n1[s])
    logit(pi1[s]) <- alpha[s]
    r2[s] ~ dbin(pi2[s], n2[s])
    logit(pi2[s]) <- alpha[s] + delta[s]
    alpha[s] ~ dnorm(mu.beta, tau.beta)   # now exchangeable with beta_h
    delta[s] ~ dnorm(mu.delta, tau.delta)
  }
  mu.delta    ~ dnorm(0, 0.1)
  sigma.delta ~ dexp(1.65); tau.delta <- pow(sigma.delta, -2)
  or          <- exp(mu.delta)
  for (h in 1:H) {
    y[h] ~ dbin(theta[h], m[h])
    logit(theta[h]) <- beta[h]
    beta[h] ~ dnorm(mu.beta, tau.beta)
  }
  mu.beta    ~ dnorm(0, 0.1)
  sigma.beta ~ dexp(1.65); tau.beta <- pow(sigma.beta, -2)
  logit(p1) <- mu.beta
  logit(p2) <- logit(p1) + mu.delta
  c.inf ~ dnorm(16.78, 0.1826)
  l     ~ dlnorm(2.075, 17.3051)
}

flu2 = jags(
  data               = dataflu,
  parameters.to.save = c("p1", "p2", "or", "c.inf", "l"),
  inits = function() {
    list(alpha = rnorm(dataflu$S, 0, .5), delta = rnorm(dataflu$S, 0, .5))
  },
  n.chains = 2, n.iter = 25000, n.burnin = 5000, n.thin = 10,
  DIC        = TRUE, pD = TRUE,
  model.file = influenza2
)

# --- Variant 2: baseline-risk adjustment for the log OR ---

#' Adds a term gamma*(alpha_s - mean(alpha)) to the treatment arm linear predictor.
#' The adjusted log OR is delta_s^* + gamma*(alpha_s - alpha-bar), where
#' delta_s^* is the unadjusted log OR.
#' mu.delta is thus the pooled effect for a study with an average baseline risk.
#' If gamma ≈ 0 (95% interval spans 0), the adjustment is immaterial.

influenza3 = function() {
  for (s in 1:S) {
    r1[s] ~ dbin(pi1[s], n1[s])
    logit(pi1[s]) <- alpha[s]
    r2[s] ~ dbin(pi2[s], n2[s])
    logit(pi2[s]) <- alpha[s] + delta[s]
    alpha[s]      ~ dnorm(0, 0.1)
    delta[s]      <- delta.unadj[s] + gamma * (alpha[s] - mean(alpha[]))
    delta.unadj[s] ~ dnorm(mu.delta, tau.delta)
  }
  gamma       ~ dnorm(0, 0.01)
  mu.delta    ~ dnorm(0, 0.1)
  sigma.delta ~ dexp(1.65); tau.delta <- pow(sigma.delta, -2)
  or          <- exp(mu.delta)
  for (h in 1:H) {
    y[h] ~ dbin(theta[h], m[h])
    logit(theta[h]) <- beta[h]
    beta[h] ~ dnorm(mu.beta, tau.beta)
  }
  mu.beta    ~ dnorm(0, 0.1)
  sigma.beta ~ dexp(1.65); tau.beta <- pow(sigma.beta, -2)
  logit(p1) <- mu.beta
  logit(p2) <- logit(p1) + mu.delta
  c.inf ~ dnorm(16.78, 0.1826)
  l     ~ dlnorm(2.075, 17.3051)
}

flu3 = jags(
  data               = dataflu,
  parameters.to.save = c("p1", "p2", "or", "c.inf", "l", "gamma"),
  inits = function() {
    list(
      alpha       = rnorm(dataflu$S, 0, .5),
      delta.unadj = rnorm(dataflu$S, 0, .5)
    )
  },
  n.chains = 2, n.iter = 25000, n.burnin = 5000, n.thin = 10,
  DIC        = TRUE, pD = TRUE,
  model.file = influenza3
)

#' Coefficient plot comparing the three model variants (p1, p2, OR).
pl1 = bmhe::coefplot(flu,  parameter = c("p1", "p2", "or"))$data
pl2 = bmhe::coefplot(flu2, parameter = c("p1", "p2", "or"))$data
pl3 = bmhe::coefplot(flu3, parameter = c("p1", "p2", "or"))$data

pl1 |>
  mutate(
    model     = "No-pooling on alpha_s",
    Parameter = case_when(
      Parameter == "p1"  ~ "p1",
      Parameter == "p2"  ~ "p2",
      TRUE               ~ "OR"
    )
  ) |>
  bind_rows(pl2 |> mutate(
    model     = "Pooling all baseline risk: alpha_s, beta_h exchangeable",
    Parameter = case_when(
      Parameter == "p1" ~ "p1", Parameter == "p2" ~ "p2", TRUE ~ "OR"
    )
  )) |>
  bind_rows(pl3 |> mutate(
    model     = "Baseline adjustment on OR",
    Parameter = case_when(
      Parameter == "p1" ~ "p1", Parameter == "p2" ~ "p2", TRUE ~ "OR"
    )
  )) |>
  ggplot(aes(mean, Parameter)) +
  geom_point(
    aes(col = model, shape = model),
    position  = position_dodge2(width = .35),
    key_glyphs = draw_key_point
  ) +
  geom_linerange(
    aes(xmin = low, xmax = upp, col = model),
    position = position_dodge2(width = .35)
  ) +
  labs(x = "Posterior distributions for the main parameters", shape = "") +
  theme(legend.position = c(.65, .45), legend.background = element_blank()) +
  labs(col = "") + ylab("") +
  scale_color_manual(values = c("#FF7F0E", "#2CA02C", "#1F77B4"))


# ==============================================================================
# SECTION: Cost-effectiveness analysis — influenza
# ==============================================================================

#' Computes population average effects and costs from the base-case flu model,
#' applying the decision tree structure in the influenza example:
#'
#' Effects (measured as negative expected illness-days):
#'   mu.e[1] = -l * p1   (status quo)
#'   mu.e[2] = -l * p2   (NIs prophylaxis)
#'
#' Costs:
#'   mu.c[1] = (1 - p1)*c.gp + p1*(c.gp + c.inf)
#'   mu.c[2] = (1 - p2)*(c.gp + c.ni) + p2*(c.gp + c.ni + c.inf)

library(BCEA)

n.sims = flu$BUGSoutput$n.sims
mu.e   = mu.c = matrix(NA, n.sims, 2)

# Fixed cost parameters (same as above)
c.gp = 19
c.ni = 2.40 * 6 * 7 * 1.175

# Extract posterior simulations
p1    = flu$BUGSoutput$sims.list$p1
p2    = flu$BUGSoutput$sims.list$p2
c.inf = flu$BUGSoutput$sims.list$c.inf
l     = flu$BUGSoutput$sims.list$l

# Population average costs
mu.c[, 1] = (1 - p1) * c.gp + p1 * (c.gp + c.inf)
mu.c[, 2] = (1 - p2) * (c.gp + c.ni) + p2 * (c.gp + c.ni + c.inf)

# Population average effects (days of influenza avoided; negative sign for direction)
mu.e[, 1] = -l * p1
mu.e[, 2] = -l * p2

#' BCEA economic analysis at WTP = £1000 per day of influenza avoided.
m = bcea(mu.e, mu.c, interventions = c("Status quo", "NIs prophylaxis"), ref = 2)
summary(m, wtp = 1000)

#' Posterior distribution of the incremental effect Delta_e = mu.e[2] - mu.e[1].
#' Displays histogram and empirical CDF to read off Pr(Delta_e > threshold).
effects = mu.e |>
  as_tibble(.name_repair = ~c("Status quo", "NIs")) |>
  mutate(Difference = NIs - `Status quo`)
costs = mu.c |>
  as_tibble(.name_repair = ~c("Status quo", "NIs")) |>
  mutate(Difference = NIs - `Status quo`)

e = ecdf(effects$Difference)
plot_eff = effects |>
  ggplot(aes(Difference)) +
  geom_bar(stat = "bin", bins = 40, col = "black", fill = "darkgrey", linewidth = .26) +
  xlim(0, 3)
maxPlotY = max(ggplot_build(plot_eff)$data[[1]]$y)

plot_eff +
  stat_ecdf(aes(y = ..y.. * maxPlotY)) +
  scale_y_continuous(
    name     = "Posterior distribution",
    sec.axis = sec_axis(transform = ~. / maxPlotY, name = "Cumulative distribution")
  ) +
  xlab("Number of days of influenza avoided") +
  geom_point(aes(x = 1, y = e(1) * maxPlotY)) +
  geom_point(aes(x = 2, y = e(2) * maxPlotY)) +
  annotate("text", x = 1, y = e(1) * maxPlotY,
           label = paste0("Pr(Difference > 1) = ", (1 - e(1)) |> round(3)),
           vjust = 1.5, hjust = -.025) +
  annotate("text", x = 2, y = e(2) * maxPlotY,
           label = paste0("Pr(Difference > 2) = ", (1 - e(2)) |> round(3)),
           vjust = 1.5, hjust = -.025) +
  geom_segment(
    aes(x = -Inf, xend = Inf, y = maxPlotY, yend = maxPlotY),
    col = "darkgrey", linetype = 2, linewidth = .26
  )
