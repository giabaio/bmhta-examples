#' ---
#' title: "Indirect treatment comparisons (Network Meta-Analysis)"
#' desc:  "Bayesian NMA of the smoking cessation dataset (S=24 studies,
#'         T=4 interventions). Covers: data wrangling to long format; prior
#'         checks for NMA parameters; fixed-effect and random-effect NMA in
#'         JAGS; consistency checks via direct vs indirect evidence; a
#'         no-pooling heterogeneity model; DIC comparison; prior visualisation
#'         for the Half-Cauchy; and a full cost-effectiveness analysis using
#'         life-years gained and BCEA."
#' ---

library(tidyverse)
library(R2jags)
library(BCEA)


# ==============================================================================
# SECTION: Data — smoking cessation (S=24 studies, T=4 interventions)
# ==============================================================================

#' The dataset is distributed as a list object (smoke.list) with elements:
#'   r  — matrix of quitter counts   (S × T_max, NA where arm not present)
#'   n  — matrix of sample sizes      (S × T_max, NA where arm not present)
#'   t  — matrix of treatment indices (S × T_max, NA where arm not present)
#'   NT — total number of treatment arms
#'
#' Interventions: 1=No intervention (A), 2=Self-help (B),
#'                3=Individual counselling (C), 4=Group counselling (D)
#'
#' Load from wherever the file is stored on your system.
load("data/nma/smoke.Rdata")   # adjust path as needed; provides smoke.list

#' Visualise the network structure (nodes = interventions, edges = studies).
ggplot() +
  annotate("text", x = 1,  y =  2, label = "A: No intervention") +
  annotate("text", x = -2, y =  1, label = "B: Self-help") +
  annotate("text", x = 4,  y = -1, label = "C: Individual counselling") +
  annotate("text", x = 0,  y = -2, label = "D: Group counselling") +
  xlim(-3, 5) + ylim(-2, 2) + theme_classic() +
  geom_segment(aes(x =  .85, y = 1.9,  xend = -2,   yend = 1.1)) +
  geom_segment(aes(x = 1.25, y = 1.9,  xend =  4,   yend = -.9)) +
  geom_segment(aes(x =  1,   y = 1.9,  xend =  0,   yend = -1.9)) +
  geom_segment(aes(x = -2,   y =  .9,  xend = -.25, yend = -1.9)) +
  geom_segment(aes(x =  4,   y = -1.1, xend =  .85, yend = -1.9)) +
  geom_segment(aes(x = -1.25, y = .92, xend =  3.0, yend = -.9)) +
  geom_label(aes(x = -.5,  y =  1.5,  label = "3 trials"),  fill = "gray80") +
  geom_label(aes(x =  2.6, y =   .5,  label = "15 trials"), fill = "gray80") +
  geom_label(aes(x =  2.6, y = -1.5,  label = "4 trials"),  fill = "gray80") +
  geom_label(aes(x = -.25, y =   .5,  label = "1 trial"),   fill = "gray80") +
  geom_label(aes(x = -1,   y = -.75,  label = "2 trials"),  fill = "gray80") +
  geom_label(aes(x =  .85, y =  1.25, label = "2 trials"),  fill = "gray80") +
  annotate("text", x = -.5,  y =  1.5, label = "N=2867",  vjust = 2.5, size = 3) +
  annotate("text", x = -1,   y = -.75, label = "N=441",   vjust = 2.5, hjust = .85, size = 3) +
  annotate("text", x = -.25, y =   .5, label = "N=255",   vjust = 2.5, hjust = .75, size = 3) +
  annotate("text", x =  2.6, y = -1.5, label = "N=764",   vjust = 2.5, size = 3) +
  annotate("text", x =  2.6, y =   .5, label = "N=12846", vjust = 2.5, hjust = .85, size = 3) +
  annotate("text", x =  .8,  y =  1.15, label = "N=318",  vjust = 1.8, hjust = -.2, size = 3) +
  theme(
    axis.line.x = element_blank(), axis.line.y = element_blank(),
    axis.ticks   = element_blank(), axis.text    = element_blank()
  ) + xlab("") + ylab("")


# ==============================================================================
# SECTION: Data wrangling — wide to long format
# ==============================================================================

#' The JAGS NMA code uses a "long-format" dataset: one row per study-arm,
#' with ancillary index variables (s, t, b) that JAGS uses as nested indices.
#' This avoids needing a fully rectangular data structure when studies have
#' different numbers of arms or different arm combinations.

colnames(smoke.list$r) = paste0("r", 1:smoke.list$NT)
colnames(smoke.list$n) = paste0("n", 1:smoke.list$NT)

smoking_data = (smoke.list$r |> as_tibble() |> mutate(s = row_number()) |>
  pivot_longer(cols = starts_with("r"), values_to = "r") |>
  drop_na()) |>
  bind_cols(
    smoke.list$n |> as_tibble() |> mutate(s = row_number()) |>
      pivot_longer(cols = starts_with("n"), values_to = "n") |>
      drop_na() |> select(-c(name, s))
  ) |>
  bind_cols(
    smoke.list$t |> as_tibble() |> mutate(s = row_number()) |>
      pivot_longer(cols = starts_with("t"), values_to = "t") |>
      drop_na() |> select(-c(name, s))
  ) |>
  group_by(s) |>
  mutate(b = min(t)) |>   # baseline = arm with the lowest treatment index
  ungroup() |>
  select(-name)

#' The 24 studies expand into 50 rows (2 arms each, plus 2 three-armed studies).
#' Key variables:
#'   s  = study index
#'   r  = quitters in this arm
#'   n  = participants in this arm
#'   t  = treatment index for this arm
#'   b  = baseline treatment index for this study

#' Shared data list for all three models.
data.list = list(
  t = smoking_data$t,
  s = smoking_data$s,
  r = smoking_data$r,
  n = smoking_data$n,
  b = smoking_data$b,
  N = nrow(smoking_data),
  T = max(smoking_data$t),
  S = max(smoking_data$s)
)


# ==============================================================================
# SECTION: Prior checks
# ==============================================================================

#' Prior predictive check: confirm that d[t] ~ Normal(0, sd=1) implies a
#' sensible range for the ORs on the natural scale.
#' Even this "tight-looking" prior allows ORs well above 7 — wide enough.
rnorm(100000, 0, 1) |> exp() |> bmhe::stats()


# ==============================================================================
# SECTION: Helper functions for direct-evidence and heterogeneity analysis
# ==============================================================================

#' direct_evidence(): fixed-effect model restricted to studies that directly
#' compare a specific pair (Active vs Baseline).  Used to check consistency
#' between direct evidence alone and the NMA-pooled indirect estimates.
direct_evidence = function() {
  for (i in 1:N) {
    r[i] ~ dbin(pi[i], n[i])
    logit(pi[i]) <- alpha[s[i]] + delta[s[i], t[i]]
    delta[s[i], t[i]] <- d[t[i]] - d[b[i]]
  }
  for (s in 1:S) {
    alpha[s] ~ dnorm(0, 0.01)
  }
  d[1] <- 0
  for (t in 2:4) {
    d[t] ~ dnorm(0, 1)
  }
  or <- exp(d[Active] - d[Baseline])   # single OR for the chosen comparison
}

#' heterogeneity(): no-pooling model — study- and treatment-specific d[s,t],
#' used to visualise raw heterogeneity across individual studies.
heterogeneity = function() {
  for (i in 1:N) {
    r[i] ~ dbin(pi[i], n[i])
    logit(pi[i]) <- alpha[s[i]] + delta[s[i], t[i]]
    delta[s[i], t[i]] <- d[s[i], t[i]] - d[s[i], b[i]]
  }
  for (s in 1:S) {
    alpha[s] ~ dnorm(0, 0.01)
    d[s, 1] <- 0
    for (j in 2:4) {
      d[s, j] ~ dnorm(0, 1)
    }
    or[s] <- exp(d[s, Active] - d[s, Baseline])
  }
}

#' make_data_no_pooling(): filters the long-format data to one pairwise
#' comparison and runs either direct_evidence() or heterogeneity().
#'
#' type = "direct":       returns pooled OR across studies for that pair
#' type = "heterogeneity": returns study-specific ORs for that pair
make_data_no_pooling = function(type = "direct") {
  comparisons = tibble(A = c(2, 3, 3, 4, 4, 4), B = c(1, 1, 2, 1, 2, 3))
  ylabs = c(
    "B: Self-help / A: None",
    "C: Individual / A: None",
    "C: Individual / B: Self-help",
    "D: Group / A: None",
    "D: Group / B: Self-help",
    "D: Group / C: Individual"
  )

  out = list()
  for (i in 1:nrow(comparisons)) {
    A = comparisons$A[i]
    B = comparisons$B[i]

    data_direct = smoking_data |>
      group_by(s) |>
      mutate(comparison = paste(t, collapse = ",")) |>
      ungroup() |>
      dplyr::filter(grepl(A, comparison), grepl(B, comparison), t %in% c(A, B)) |>
      select(r, n, t, b, s) |>
      group_by(s) |> mutate(ss = cur_group_id()) |>
      ungroup() |> select(-s) |> rename(s = ss) |>
      as.list()
    data_direct$N        = length(data_direct$r)
    data_direct$S        = max(data_direct$s)
    data_direct$Baseline = B
    data_direct$Active   = A

    model_fn = if (type == "direct") direct_evidence else heterogeneity
    m = jags(
      data = data_direct, parameters.to.save = c("or"), n.thin = 4,
      inits = NULL, n.chains = 2, n.burnin = 1000, n.iter = 10000,
      model.file = model_fn
    )

    out[[i]] = m$BUGSoutput$summary[, c("mean", "sd", "2.5%", "97.5%")] |>
      as_tibble() |>
      mutate(Parameter = rownames(m$BUGSoutput$summary[, c("mean","sd","2.5%","97.5%")])) |>
      dplyr::filter(Parameter != "deviance")

    if (type == "heterogeneity") {
      out[[i]] = out[[i]] |>
        mutate(
          comparison = ylabs[i],
          N = data_direct |> as_tibble() |> select(n, s) |>
            group_by(s) |> summarise(N = sum(n)) |> pull(N),
          s = smoking_data |> group_by(s) |>
            mutate(comparison = paste(t, collapse = ",")) |> ungroup() |>
            dplyr::filter(
              grepl(A, comparison), grepl(B, comparison), t %in% c(A, B)
            ) |> select(s) |> unique() |> pull(s)
        )
    }
  }

  if (type == "direct") {
    out = out |> bind_rows() |>
      mutate(
        comparison = ylabs,
        Parameter  = paste0("or[", comparisons$A, ",", comparisons$B, "]"),
        num        = row_number()
      ) |>
      select(Parameter, mean, sd, `2.5%`, `97.5%`, num)
  } else {
    out = out |> bind_rows()
  }
  return(out)
}


# ==============================================================================
# SECTION: Model 1 — Fixed-effect NMA (complete pooling)
# ==============================================================================

#' Assumes a single common log OR d[t] for each treatment vs the baseline.
#' delta[s,t] = d[t] - d[b[s]] is the study-arm-specific incremental effect,
#' which in the fixed-effect model is purely deterministic given d[t] and d[b[s]].
#'
#' The absolute probability of quitting under each intervention is derived
#' using an informative prior for the baseline log-odds of quitting:
#'   rho ~ Normal(-2.6, precision=6.925)   [implies sd ≈ 0.38]
#' This encodes: P(quit | no intervention) ≈ 7%, upper limit ≈ 14%.
#'
#' All pairwise ORs are monitored in both directions:
#'   or[A,B] = exp(d[A] - d[B]);   or[B,A] = 1/or[A,B]

smoking_fe = function() {
  for (i in 1:N) {
    r[i] ~ dbin(pi[i], n[i])
    logit(pi[i]) <- alpha[s[i]] + delta[s[i], t[i]]
    delta[s[i], t[i]] <- d[t[i]] - d[b[i]]
  }
  for (s in 1:S) {
    alpha[s] ~ dnorm(0, 0.01)   # study baseline, Normal(0, sd=10) on logit scale
  }
  d[1] <- 0                      # reference arm has no incremental effect
  for (t in 2:T) {
    d[t] ~ dnorm(0, 1)           # log OR, Normal(0, sd=1): wide on OR scale
  }
  # All pairwise ORs
  for (A in 1:(T - 1)) {
    for (B in (A + 1):T) {
      or[A, B] <- exp(d[A] - d[B])
      or[B, A] <- 1 / or[A, B]
    }
  }
  # Baseline probability of quitting with no intervention (external evidence)
  rho ~ dnorm(-2.6, 6.925208)    # precision = 1/0.38^2 ≈ 6.925
  # Absolute probability of quitting under each intervention
  for (t in 1:T) {
    logit(p[t]) <- rho + d[t]
  }
}

m_fe = jags(
  data               = data.list,
  parameters.to.save = c("d", "or", "p"),
  inits = NULL, n.chains = 2, n.burnin = 5000, n.iter = 25000, n.thin = 10,
  model.file = smoking_fe
)


# ==============================================================================
# SECTION: Direct vs indirect evidence comparison
# ==============================================================================

#' Runs make_data_no_pooling("direct") to get OR estimates from studies that
#' directly compare each pair, then overlays them with the fixed-effect NMA
#' (which uses all evidence, including indirect).
#'
#' Comparisons with large sample sizes (e.g. C vs A) are very stable; the
#' discrepancy is most visible for pairs with sparse direct data (e.g. C vs B).

ylabs = c(
  "B: Self-help / A: None (N=2867)",
  "C: Individual / A: None (N=12846)",
  "C: Individual / B: Self-help (N=255)",
  "D: Group / A: None (N=318)",
  "D: Group / B: Self-help (N=441)",
  "D: Group / C: Individual (N=764)"
)

select = c("or[2,1]", "or[3,1]", "or[3,2]", "or[4,1]", "or[4,2]", "or[4,3]")
toplot = bmhe::coefplot(m_fe, parameter = select)$data

#' Direct-evidence only estimates (no pooling across studies for each pair).
out1 = make_data_no_pooling(type = "direct") |> rename(low = `2.5%`, upp = `97.5%`)

toplot |> mutate(model = "Fixed effect (indirect)", Parameter = ylabs) |>
  bind_rows(out1 |> mutate(model = "No pooling effect (direct)", Parameter = ylabs)) |>
  ggplot(aes(mean, Parameter)) +
  geom_point(aes(col = model, shape = model),
             position = position_dodge2(width = .25)) +
  geom_linerange(aes(xmin = low, xmax = upp, col = model),
                 position = position_dodge2(width = .25)) +
  geom_vline(xintercept = 1, col = "gray50") +
  labs(x = "Odds Ratio", shape = "") +
  theme(legend.position = c(.75, .1), legend.background = element_blank()) +
  labs(col = "") + ylab("") + xlim(0, 6) +
  scale_color_manual(values = c("#FF7F0E", "#000000"))


# ==============================================================================
# SECTION: Heterogeneity across individual studies
# ==============================================================================

#' Runs the no-pooling heterogeneity model separately for each pairwise
#' comparison, producing study-specific OR estimates.
#' The resulting forest plot (log scale x-axis, faceted by comparison) reveals
#' substantial between-study heterogeneity, motivating the random-effects model.

out2 = make_data_no_pooling(type = "heterogeneity") |> rename(low = `2.5%`, upp = `97.5%`)

out2 |>
  mutate(Parameter = as.factor(paste0("s=", s, " N=", N))) |>
  ggplot(aes(mean, fct_reorder(Parameter, s, .desc = TRUE))) +
  geom_linerange(aes(xmin = low, xmax = upp)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_point() +
  labs(x = "Odds Ratio (log scale axis)") +
  facet_grid(comparison ~ ., scales = "free", space = "free") +
  ylab("") +
  theme(strip.text.y = element_text(angle = 0)) +
  scale_x_continuous(
    trans   = "log",
    labels  = c(.2, .5, 1, 2, 5, 10, 25),
    breaks  = c(.2, .5, 1, 2, 5, 10, 25)
  )


# ==============================================================================
# SECTION: Prior for the Half-Cauchy random-effects SD
# ==============================================================================

#' The Half-Cauchy(0, lambda) prior for sigma is equivalent to a truncated
#' Student-t distribution with df=1, mu=0 and precision = lambda^{-2}, 
#' restricted to sigma > 0 (via the T(0,) truncation operator in JAGS).
#'
#' We use Half-Cauchy(0, 2) -> precision = 0.25 in the JAGS dt() call.
#' A helper function ttrunc() computes the density of the mean-precision t
#' in R (using base R's dt() and pt(), which use the standardised parameterisation).

ttrunc = function(x, mu = 0, tau = 1, df = 1, left = -Inf, right = Inf) {
  f = sqrt(tau) * dt(((x - mu) * sqrt(tau)), df = df) /
    (pt(((right - mu) * sqrt(tau)), df = df) - pt(((left - mu) * sqrt(tau)), df = df))
  return(f)
}

#' Pr(sigma > 6) under the Half-Cauchy(0,2) prior (forward sampling check).
p6 = 1 - pt((6 * sqrt(.25)), df = 1) |> round(3)
cat("Pr(sigma > 6) under Half-Cauchy(0,2):", p6, "\n")

#' Visualise Half-Cauchy(0,2) vs Exponential(0.4) — the two priors are
#' closely aligned, confirming the PC-prior interpretation of the Half-Cauchy.
tibble(x = seq(0, 20, .01), y = ttrunc(seq(0, 20, .01), mu = 0, tau = 0.25, left = .0001)) |>
  ggplot(aes(x, y, color = "hc", linetype = "solid")) + geom_line() +
  stat_function(fun = dexp, args = list(.4), aes(colour = "exp", linetype = "dashed")) +
  xlim(0, 12) +
  xlab("sigma") + ylab("") +
  theme(legend.position = c(.65, .75)) +
  scale_color_manual(
    name   = "",
    values = c("hc" = "#E69F00", "exp" = "#0072B2"),
    labels = c("Exponential(0.4)", "Half Cauchy(0,2)")
  ) +
  theme(legend.background = element_blank()) +
  scale_linetype_manual(
    values = c("solid", "dashed"), name = "",
    labels = c("Exponential(0.4)", "Half Cauchy(0,2)")
  )

#' Note on R vs JAGS parameterisation of the t distribution:
#' In R, dt(y, df) assumes mu=0 and tau=1.
#' In JAGS, dt(mu, tau, df) uses the mean-precision parameterisation.
#' If y ~ dt(mu, tau, df) in JAGS, then z = (y - mu)*sqrt(tau) ~ dt(0, 1, df) in R.
#' Density and CDF in R:
#'   density  <- sqrt(tau) * dt(sqrt(tau) * (y - mu), df)
#'   cum_dist <- pt(sqrt(tau) * (y - mu), df)


# ==============================================================================
# SECTION: Model 2 — Random-effects NMA (partial pooling)
# ==============================================================================

#' The random-effects NMA allows study-specific log ORs delta[s,t] to vary
#' around the pooled estimate mu[s,t] = d[t] - d[b[s]], with a common
#' standard deviation sigma across all treatment comparisons.
#'
#' Key structural addition vs fixed-effect model:
#'   delta[s,t] ~ Normal(mu[s,t], tau)   [stochastic, not deterministic]
#'   sigma ~ Half-Cauchy(0, 2)            [via JAGS truncated dt]
#'   tau   = sigma^{-2}
#'
#' IMPORTANT: the truncation operator T(0,) cannot appear directly inside an
#' R function because it breaks R syntax.  The workaround is to write %_%T(0,)
#' which R2jags/R2WinBUGS parses out to T(0,) in the JAGS model file.

smoking_re = function() {
  for (i in 1:N) {
    r[i] ~ dbin(pi[i], n[i])
    logit(pi[i]) <- alpha[s[i]] + delta[s[i], t[i]]
    delta[s[i], t[i]] ~ dnorm(mu[s[i], t[i]], tau)   # stochastic RE
    mu[s[i], t[i]]    <- d[t[i]] - d[b[i]]
  }
  for (s in 1:S) {
    alpha[s] ~ dnorm(0, 0.01)
  }
  d[1] <- 0
  for (t in 2:T) {
    d[t] ~ dnorm(0, 1)
  }
  # Half-Cauchy(0,2) prior for the common RE standard deviation
  # In JAGS: dt(mu, tau, df) with tau=precision; precision=0.25 => scale lambda=2
  sigma ~ dt(0, 0.25, 1) %_% T(0, )
  tau   <- pow(sigma, -2)

  # All pairwise ORs
  for (A in 1:(T - 1)) {
    for (B in (A + 1):T) {
      or[A, B] <- exp(d[A] - d[B])
      or[B, A] <- 1 / or[A, B]
    }
  }
  # Baseline probability of quitting with no intervention
  rho ~ dnorm(-2.6, 6.925208)
  for (t in 1:T) {
    logit(p[t]) <- rho + d[t]
  }
}

m_re = jags(
  data               = data.list,
  parameters.to.save = c("d", "or", "p", "sigma"),
  inits = NULL, n.chains = 2, n.burnin = 5000, n.iter = 25000, n.thin = 10,
  model.file = smoking_re
)

# --- Diagnostics ---
bmhe::diagplot(m_re)
bmhe::diagplot(m_re, what = "n.eff", label = TRUE) + ylim(0, 4250)


# ==============================================================================
# SECTION: DIC comparison and effective parameter count
# ==============================================================================

#' Computes pD for the random-effects model and derives the % shrinkage
#' relative to the nominal number of stochastic nodes.
#' This measures how much correlation / borrowing the RE structure induces.

pD = rjags::dic.samples(m_re$model, n.iter = 1000)$penalty |> sum()

#' count_nodes() extracts the stochastic node count from a JAGS object.
#' If not available in your version of bmhe, the nominal count can be read
#' directly from the JAGS compile output (Unobserved stochastic nodes: ...).
#' For this model the nominal count is ~51 (16 alpha + 34 delta + 1 sigma).
cat("Fixed effect DIC: ",  m_fe$BUGSoutput$DIC |> round(2), "\n")
cat("Random effect DIC:", m_re$BUGSoutput$DIC |> round(2), "\n")
cat("pD (RE model):", pD |> round(2), "\n")


# ==============================================================================
# SECTION: Full comparison — direct, fixed effect, random effect
# ==============================================================================

#' Overlays all three model types on the same forest plot.
#' The random-effects model generally produces wider intervals (proper
#' accounting of heterogeneity) compared to the artificially tight fixed-effect
#' intervals, and acts as a compromise between direct and fixed-effect estimates.

toplot_fe = bmhe::coefplot(m_fe, parameter = select)$data
toplot_re = bmhe::coefplot(m_re, parameter = select)$data

toplot_fe |> mutate(model = "Fixed effect (indirect)", Parameter = ylabs) |>
  bind_rows(out1 |> mutate(model = "No pooling (direct)", Parameter = ylabs)) |>
  bind_rows(toplot_re |> mutate(model = "Random effects", Parameter = ylabs)) |>
  ggplot(aes(mean, Parameter)) +
  geom_point(aes(col = model, shape = model),
             position = position_dodge2(width = .35)) +
  geom_linerange(aes(xmin = low, xmax = upp, col = model),
                 position = position_dodge2(width = .35)) +
  geom_vline(xintercept = 1, col = "gray50") +
  labs(x = "Odds Ratio", shape = "") +
  theme(legend.position = c(.75, .12), legend.background = element_blank()) +
  labs(col = "") + ylab("") + xlim(0, 6) +
  scale_color_manual(values = c("#FF7F0E", "#000000", "#1F77B4"))


# ==============================================================================
# SECTION: Absolute probability of quitting smoking
# ==============================================================================

#' Compares posterior distributions for p[t] (absolute probability of quitting)
#' between the fixed and random effects models, labelled by intervention.
#' The random-effects model gives wider intervals, reflecting heterogeneity.
#' The "No intervention" arm is essentially identical in both models because
#' it is anchored to the informative prior on rho.

toplot_fe = bmhe::coefplot(m_fe, parameter = "p")$data
toplot_re = bmhe::coefplot(m_re, parameter = "p")$data

int_labels = c(
  "p[1]" = "No intervention", "p[2]" = "Self-help",
  "p[3]" = "Individual counselling", "p[4]" = "Group counselling"
)

toplot = toplot_fe |> mutate(model = "Fixed effects") |>
  bind_rows(toplot_re |> mutate(model = "Random effects")) |>
  mutate(
    Parameter = str_replace_all(Parameter, "p\\[1\\]", "No intervention"),
    Parameter = str_replace_all(Parameter, "p\\[2\\]", "Self-help"),
    Parameter = str_replace_all(Parameter, "p\\[3\\]", "Individual counselling"),
    Parameter = str_replace_all(Parameter, "p\\[4\\]", "Group counselling"),
    Parameter = factor(Parameter, levels = c(
      "No intervention", "Self-help",
      "Individual counselling", "Group counselling"
    ))
  )

toplot |> ggplot(aes(mean, Parameter)) +
  geom_point(aes(col = model, shape = model),
             position = position_dodge2(width = .25)) +
  geom_linerange(aes(xmin = low, xmax = upp, col = model),
                 position = position_dodge2(width = .25)) +
  labs(x = "Probability of quitting smoking") +
  theme(legend.position = c(.75, .15), legend.background = element_blank()) +
  labs(col = "", shape = "") + ylab("") +
  scale_color_manual(values = c("#FF7F0E", "#1F77B4")) +
  scale_shape_manual(values = c(16, 15))


# ==============================================================================
# SECTION: Cost-effectiveness analysis
# ==============================================================================

#' Converts absolute quitting probabilities to life-years gained (LYG), using
#' external evidence from Mamun et al. (2004) via the Framingham Heart Study:
#'   Men:   LYG ~ Normal(8.66, sd=0.495)  [sd from 95% CI = (7.61, 9.63)]
#'   Women: LYG ~ Normal(7.59, sd=0.679)  [sd from 95% CI = (6.33, 8.92)]
#'
#' Gender mix based on ONS 2024 UK smoking statistics:
#'   3.4M men : 2.6M women -> weights 0.567 and 0.433.
#'
#' Effect for each intervention: e[,t] = p[t] * l
#' where l = 0.567*l_m + 0.433*l_f is the weighted average LYG.

lm = rnorm(m_re$BUGSoutput$n.sims, 8.66, 0.495)  # LYG for men
lf = rnorm(m_re$BUGSoutput$n.sims, 7.59, 0.679)  # LYG for women
l  = 0.567 * lm + 0.433 * lf                      # weighted average LYG

#' e is a matrix: n.sims rows × 4 columns (one per intervention).
e = m_re$BUGSoutput$sims.list$p * l

labs = c("No intervention", "Self-help", "Individual counselling", "Group counselling")

#' Summary of the intervention-specific posterior LYG distributions.
e |> bmhe::stats() |> (\(x) { rownames(x) = labs; x })()

# --- Cost parameters ---

#' Cost point estimates (£) from BCEA (Table 4.2).
#' Items: NRT = 35 nicotine patches × £1.30 = £45.50
#'   No intervention:        £0
#'   Self-help (B):          £45.50 (NRT only)
#'   Individual counselling: £45.50 + 5 clinic visits × £10 = £95.50
#'   Group counselling:      £45.50 + 5 clinic visits × £19.46 = £142.80
#'
#' Uncertainty is modelled as Gamma(shape, rate) with CV = sd/mean = 0.3.
#' bmhe::gammaPar(mean, sd) returns the shape and rate for a Gamma matching
#' the given natural-scale mean and sd.

bmhe::gammaPar(45.5,  0.3 * 45.5)   # Self-help
bmhe::gammaPar(95.5,  0.3 * 95.5)   # Individual counselling
bmhe::gammaPar(142.8, 0.3 * 142.8)  # Group counselling

#' Quick Monte Carlo check: verify the Gamma parameters recover the intended
#' mean and sd for self-help.
bmhe::gammaPar(45.5, .3 * 45.5) |>
  (\(x) rgamma(100000, shape = x$shape, rate = x$rate))() |>
  bmhe::stats()

#' Build the cost simulation matrix (n.sims × 4).
#' No intervention has zero cost with no variance.
c = matrix(
  c(
    rep(0, m_re$BUGSoutput$n.sims),                                         # A: £0
    rgamma(m_re$BUGSoutput$n.sims, shape = 11.111, rate = 0.244),           # B
    rgamma(m_re$BUGSoutput$n.sims, shape = 11.111, rate = 0.116),           # C
    rgamma(m_re$BUGSoutput$n.sims, shape = 11.111, rate = 0.078)            # D
  ),
  nrow = m_re$BUGSoutput$n.sims, ncol = 4
)

#' BCEA economic analysis with WTP up to £4000/LYG.
#' ref=4 sets Group counselling (D) as the reference intervention.
m_cea = bcea(e, c, ref = 4, Kmax = 4000, interventions = labs)
summary(m_cea)

#' Multi-intervention CEA: probability each intervention is most cost-effective
#' across the WTP range, and the cost-effectiveness efficiency frontier.
mce = multi.ce(m_cea)

cc = ceac.plot(mce, graph = "gg")
cc$data |>
  mutate(label = labs[comparison] |> factor(levels = labs)) |>
  ggplot(aes(k, ceac, col = label)) +
  geom_line(linewidth = .95, key_glyph = draw_key_label) +
  labs(title = "") +
  theme(
    legend.position.inside = c(.55, .80), legend.position = "inside",
    legend.background       = element_blank(),
    legend.text             = element_text(size = 12)
  ) +
  labs(color = "") +
  xlab("Willingness to pay") +
  ylab("Probability of most cost-effective") +
  scale_x_continuous(label = scales::comma, limits = c(0, 4050)) +
  directlabels::geom_dl(
    aes(label = comparison),
    method = list("last.points", cex = 1.2, colour = "black")
  ) +
  scale_color_manual(values = c("#000000", "#1F77B4", "#FF7F0E", "#2CA02C")) +
  guides(colour = guide_legend(override.aes = list(label = c("1", "2", "3", "4"))))

ceef.plot(mce, graph = "gg", print.summary = FALSE) +
  labs(title = "") +
  theme(
    legend.text             = element_text(size = 12),
    legend.position.inside  = c(.8, .8),
    legend.position         = "inside"
  )
