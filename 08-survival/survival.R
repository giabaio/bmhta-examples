#' ---
#' title: "Survival analysis in HTA"
#' desc:  "Bayesian survival modelling for HTA, using the digitised NICE TA174
#'         data (rituximab RFC vs FC for chronic lymphocytic leukaemia) as a
#'         running example. Covers: censoring visualisation; Kaplan-Meier
#'         estimation; parametric Bayesian survival models (Weibull, Gompertz,
#'         log-Normal, Generalised F) via survHE/rstan; hazard and cumulative
#'         hazard plots; probabilistic sensitivity analysis; model selection via
#'         DIC; separate modelling by arm; extrapolation; M-spline models via
#'         survextrap; and incorporation of external aggregated data."
#' ---

library(tidyverse)
library(survHE)
library(survextrap)


# ==============================================================================
# SECTION: Data — NICE TA174 (digitised)
# ==============================================================================

#' The TA174 dataset is bundled with survHE.  It is a digitised and
#' post-processed version of the CLL-8 trial (Hallek et al. 2010), comparing
#' RFC (rituximab + fludarabine + cyclophosphamide) against FC alone for
#' first-line treatment of chronic lymphocytic leukaemia.
#'
#' Variables:
#'   patid   — patient ID
#'   prog    — progression indicator (1=progressed, 0=censored)
#'   death   — death indicator (1=died, 0=censored)
#'   prog_t  — time to progression (or censoring), months
#'   death_t — time to death (or censoring), months
#'   treat   — treatment arm (0=FC, 1=RFC)
data(TA174, package = "survHE")

ta174 = ta174 |>
  mutate(
    treatment = case_when(treat == 0 ~ "FC", TRUE ~ "RFC") |> as.factor()
  )

#' We focus on progression as the endpoint.  Individuals who die without
#' progression, or who are censored without progressing, are treated as
#' censored in this analysis (status = 0).
data = ta174 |>
  mutate(
    id        = patid,
    from      = 1, to = 2, trans = 1,
    Tstart    = 0,
    Tstop     = prog_t,
    time      = Tstop - Tstart,
    status    = case_when(prog == 0 ~ 0, TRUE ~ 1),
    treat     = treat,
    treatment = treatment
  ) |>
  select(id, from, to, trans, Tstart, Tstop, time, status, treatment, treat)


# ==============================================================================
# SECTION: Censoring illustration
# ==============================================================================

#' Draws a "swim-lane" diagram showing five fictitious individuals, some of
#' whom experience the event (filled circle) and some of whom are censored
#' (dashed line ending before the study close).  Columns to the right show
#' the observed time t_i, event indicator d_i and whether the true event time
#' is known.

dat.cens = tibble(
  start = c(1, 4, 2, 3, 1),
  end   = c(6, 5, 4, 6.5, 5),
  obs   = c(1, 1, 0, 0, 1),
  y     = c(1, 2, 3, 4, 5)
) |>
  mutate(
    tt = case_when(
      obs == 1 ~ glue::glue("{end}"),
      TRUE     ~ glue::glue("? (>= {end})")
    )
  )

dat.cens |>
  ggplot() +
  geom_segment(aes(x = start, xend = end, y = y, yend = y)) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),  axis.line.x = element_blank()
  ) +
  geom_point(data = data.frame(x = c(6, 5, 5), y = c(1, 2, 5)),
             aes(x, y), size = 3) +
  geom_point(data = data.frame(x = c(1, 4, 2, 3, 1), y = c(1, 2, 3, 4, 5)),
             aes(x, y), size = 3, shape = 15) +
  geom_vline(xintercept = 6.5) +
  annotate("text", x = 7,    y = 6, label = "t_i",       size = 7) +
  annotate("text", x = 7.75, y = 6, label = "d_i",       size = 7) +
  annotate("text", x = 8.95, y = 6, label = "True time", size = 6) +
  annotate("text", x = 7,    y = 1:5, label = dat.cens$end, size = 5) +
  annotate("text", x = 7.75, y = 1:5, label = dat.cens$obs, size = 5) +
  annotate("text", x = 8.95, y = 1:5, label = dat.cens$tt,  size = 5) +
  xlab("follow-up period") + ylab("Individuals") + xlim(0, 10)


# ==============================================================================
# SECTION: Kaplan-Meier estimation (TA174)
# ==============================================================================

#' fit.models() with a single distribution and no method= argument defaults to
#' MLE via flexsurv.  Here we use it primarily to obtain the KM curves,
#' which are stored in m$misc$km and can be plotted via plot(m, add.km=TRUE).
m_km = fit.models(Surv(time, status) ~ treatment, data, distr = "wei")

#' Plot the KM curves for the two treatment arms.
#' We remove the parametric model layer (layer 1) to show only the KM.
p_km = plot(m_km, add.km = TRUE, lab.profile = c("RFC", "FC"))
p_km$layers[[1]] = NULL
p_km +
  guides(color = "none") +
  annotate("text", x = 20, y = .73, label = "FC") +
  annotate("text", x = 30, y = .77, label = "RFC")

#' Extract the underlying KM numbers (events, at risk, censored) per arm
#' for tabulation.
kmdata = tibble(
  time     = m_km$misc$km$time,
  n_event  = m_km$misc$km$n.event,
  n_risk   = m_km$misc$km$n.risk,
  n_censor = m_km$misc$km$n.censor
) |>
  mutate(
    treatment = c(
      rep("RFC", m_km$misc$km$strata[1]),
      rep("FC",  m_km$misc$km$strata[2])
    )
  )


# ==============================================================================
# SECTION: Generalised F model — Bayesian advantage illustration
# ==============================================================================

#' The Generalised F has 3 ancillary parameters (sigma, p, q).
#' Under MLE with limited/censored data, some parameters can become
#' non-identifiable (e.g. p's 95% CI spans essentially [0, Inf]).
#' A Bayesian run with regularising priors fixes this.

#' --- MLE run (for comparison) ---
load("data/survival/data_surv.Rdata")
dat$TIME=dat$time
dat$EVENT=dat$censored
dat$treatment=factor(dat$arm,labels=c("Comparator","Intervention"))

m0 = survHE::fit.models(
  Surv(time, censored) ~ treatment,
  distr = "genf", data = dat   # 'dat' is a different dataset used in the chapter
)

#' --- Bayesian run (HMC via rstan) ---
#' Default priors in survHE: sigma ~ Gamma(0.1,0.1), log(p) ~ Normal(0,0.5),
#' q ~ Normal(0,2.5).  These regularise the otherwise-unidentifiable p.
m1 = survHE::fit.models(
  Surv(time, censored) ~ treatment,
  distr = "genf", data = dat,
  method = "hmc"
)

#' Compare: MLE returns nearly infinite CI for p; Bayesian gives finite estimate.
print(m0)
print(m1)


# ==============================================================================
# SECTION: Fitting Bayesian parametric models via survHE (TA174)
# ==============================================================================

#' Fits Weibull (AFT), Gompertz and log-Normal models simultaneously using
#' Hamiltonian Monte Carlo (Stan/rstan) via survHE.
#'
#' Default priors used by survHE:
#'   Weibull/Gompertz: beta ~ Normal(0, 5),  alpha ~ Gamma(0.1, 0.1)
#'   log-Normal:       beta ~ Normal(0, 100), log(alpha) ~ Uniform(0, 5)
#'
#' fit.models() returns a list with:
#'   models         — rstan output for each distribution
#'   model.fitting  — AIC, BIC and DIC
#'   method         — "hmc"
#'   misc           — data, formula, running time, KM object

m = fit.models(
  Surv(time, status) ~ treatment,
  data   = data,
  distr  = c("wei", "gom", "lno"),
  method = "hmc"
)

#' Print posterior summaries for each model.
print(m, mod = 1)   # Weibull
print(m, mod = 2)   # Gompertz
print(m, mod = 3)   # log-Normal

#' Print with original Stan parameter names, prior information and
#' Stan convergence metrics (Rhat, n_eff).
print(m, mod = 2, original = TRUE, digits = 2, print_priors = TRUE)


# ==============================================================================
# SECTION: Convergence diagnostics (rstan interface)
# ==============================================================================

#' Because m$models[[k]] is a rstan object, all rstan diagnostics apply.
#' traceplot() and stan_ac() are the two most useful for quick inspection.

# Traceplot for Weibull model parameters
rstan::traceplot(m$models[[1]])

# Autocorrelation plots for Weibull model parameters
rstan::stan_ac(m$models[[1]])


# ==============================================================================
# SECTION: Survival curve plots
# ==============================================================================

#' plot() in survHE returns a ggplot object — it can be further customised
#' using standard ggplot2 facilities.

#' All three models with KM overlay (observation period only).
plot(
  m,
  add.km      = TRUE,
  lab.profile = c("RFC", "FC")
) +
  theme(legend.position.inside = c(.2, .2)) +
  scale_color_manual(values = c("#000000", "#1F77B4", "#FF7F0E"))

#' Hazard functions.
plot(m, what = "hazard", lab.profile = c("RFC", "FC")) +
  theme(legend.position.inside = c(.2, .8)) +
  scale_color_manual(values = c("#000000", "#1F77B4", "#FF7F0E"))

#' Cumulative hazard functions (= -log S(t), by eq-cumhaz).
plot(m, what = "cumhazard", lab.profile = c("RFC", "FC")) +
  theme(legend.position.inside = c(.2, .8)) +
  scale_color_manual(values = c("#000000", "#1F77B4", "#FF7F0E"))


# ==============================================================================
# SECTION: Probabilistic sensitivity analysis (posterior survival bands)
# ==============================================================================

#' Setting nsim > 1 propagates the joint posterior uncertainty in the
#' parameters to a distribution of survival curves, shown as a shaded band.
#' mods=2 selects only the Gompertz model to avoid cluttering the plot.

p = plot(
  m,
  mods        = 2,
  add.km      = TRUE,
  lab.profile = c("RFC", "FC"),
  nsim        = 1000
)

#' The fourth layer is the CI ribbon around the KM; remove it to keep only
#' the KM step curve, which makes the plot more readable.
p$layers[[4]] = NULL
p + theme(legend.position.inside = c(.2, .6))


# ==============================================================================
# SECTION: Model selection (DIC)
# ==============================================================================

#' model.fit.plot() visualises DIC across all fitted models as a bar chart.
#' type="DIC" shows absolute values; scale="relative" shows % increase over
#' the best-fitting model.

model.fit.plot(m, type = "DIC")
model.fit.plot(m, type = "DIC", scale = "relative")


# ==============================================================================
# SECTION: Separate modelling by treatment arm
# ==============================================================================

#' When no single distributional assumption fits both arms equally well,
#' we can fit separate models per arm using only an intercept (~1).
#' This relaxes the PH assumption: location AND ancillary parameters are
#' allowed to differ between arms.

m.ctl = fit.models(
  Surv(time, status) ~ 1,
  data   = data |> dplyr::filter(treatment == "FC"),
  distr  = c("wei", "gom", "lno"),
  method = "hmc"
)

m.trt = fit.models(
  Surv(time, status) ~ 1,
  data   = data |> dplyr::filter(treatment == "RFC"),
  distr  = c("wei", "gom", "lno"),
  method = "hmc"
)

#' Stacked DIC comparison across both arm-specific objects.
#' stacked=TRUE groups bars by model for easier cross-arm comparison.
model.fit.plot(FC = m.ctl, RFC = m.trt, type = "DIC", stacked = TRUE)


# ==============================================================================
# SECTION: Extrapolation
# ==============================================================================

#' Extrapolate to t=180 months (15 years) by passing t=seq(0,180) to plot().
plot(
  m,
  add.km      = TRUE,
  lab.profile = c("RFC", "FC"),
  t           = seq(0, 180)
) +
  scale_color_manual(values = c("#000000", "#1F77B4", "#FF7F0E"))

#' Select best per-arm model: Weibull (FC arm) + Gompertz (RFC arm).
#' mods=c(1,5): 1st model in m.ctl (Weibull for FC) and 5th overall =
#' 2nd model in m.trt (Gompertz for RFC).
plot(
  FC = m.ctl, RFC = m.trt,
  mods   = c(1, 5),
  add.km = TRUE,
  t      = seq(0, 180)
) +
  scale_color_manual(values = c("#000000", "#1F77B4"))


# ==============================================================================
# SECTION: M-spline model via survextrap
# ==============================================================================

#' M-splines provide a flexible, non-parametric-style hazard model that avoids
#' committing to a specific distributional shape.  The hazard is expressed as
#' a weighted sum of basis functions; weights are partially pooled (penalised)
#' to prevent overfitting.
#'
#' Key survextrap workflow:
#'   1. mspline_spec(): defines the spline structure (knots, df, extrapolation)
#'   2. survextrap():   fits the Bayesian model via Stan
#'   3. survival():     extracts posterior summaries of S(t)

#' --- Step 1: spline specification ---
#' df=6 basis functions; add_knots=180 ensures extrapolation up to 15 years.
#' Knot locations are auto-selected from quantiles of observed event times.
mspline = mspline_spec(
  formula   = Surv(time, status) ~ 1,
  data      = data,
  df        = 6,
  add_knots = 180
)

#' --- Step 2: model fit ---
#' treatment is included as a covariate in the formula; under the default
#' PH assumption, a single HR is estimated (loghr in the output).
msp = survextrap(
  formula = Surv(time, status) ~ treatment,
  data    = data,
  mspline = mspline,
  chains  = 2,
  iter    = 4000
)

#' Print model summary: shows prior setup, parameter posteriors (alpha=log phi,
#' coefs=omega_k, loghr=log HR for treatment, hr=HR, hsd=smoothing SD) and
#' Stan convergence metrics.
msp

#' --- Step 3: extract survival curves ---
#' survival() returns a tidy tibble with t, treatment, mean and quantiles.
#' Using mean (not median) for consistency with survHE conventions.
S = survextrap::survival(
  msp,
  summ_fns = list(mean = mean, ~quantile(.x, probs = c(0.025, 0.975)))
)

#' --- Combine survHE and survextrap plots ---
#' Overlays Weibull + Gompertz (from survHE) with M-spline (from survextrap).
#' The treatment strata labels are recoded to match survHE's internal naming.
plot(
  m, mods = c(1, 2), add.km = TRUE, t = seq(0, 180),
  lab.profile = c("RFC", "FC")
) +
  geom_line(
    data = S |>
      mutate(
        model  = "M-spline",
        strata = case_when(
          treatment == "FC" ~ "treatmentRFC=0",
          TRUE              ~ "treatmentRFC=1"
        ) |> factor(levels = c("treatmentRFC=1", "treatmentRFC=0"))
      ),
    aes(t, mean, col = model, linetype = strata),
    linewidth = 0.9
  ) +
  scale_color_manual(values = c("#000000", "#1F77B4", "#FF7F0E"))

#' Uncertainty band for the M-spline model alone (PSA-style visualisation).
#' Wide ribbons in the extrapolation region are expected: the spline is
#' agnostic about long-run behaviour, unlike parametric models which are
#' forced by their mathematical form to eventually reach 0.
fig_mspline = plot(m, mods = c(1, 2), add.km = TRUE, t = seq(0, 180),
                   lab.profile = c("RFC", "FC")) +
  geom_line(
    data = S |>
      mutate(
        model  = "M-spline",
        strata = case_when(
          treatment == "FC" ~ "treatmentRFC=0", TRUE ~ "treatmentRFC=1"
        ) |> factor(levels = c("treatmentRFC=1", "treatmentRFC=0"))
      ),
    aes(t, mean, col = model, linetype = strata), linewidth = 0.9
  ) +
  scale_color_manual(values = c("#000000", "#1F77B4", "#FF7F0E"))

#' Extract the M-spline layer from the combined plot for a standalone PSA view.
fig_mspline$layers[[4]]$data |>
  ggplot(aes(t, mean, col = model)) +
  geom_ribbon(
    aes(ymin = `2.5%`, ymax = `97.5%`, group = strata),
    alpha = 0.2, linetype = 0, col = NA, show.legend = FALSE
  ) +
  geom_line(aes(linetype = strata, color = model), linewidth = 0.9) +
  theme_survHE() +
  theme(legend.position = c(0.6, 0.95), legend.justification = c("left", "top")) +
  labs(color = "Model", linetype = "Profile", x = "Time", y = "Survival") +
  scale_linetype_manual(
    values = c("treatmentRFC=0" = "dashed", "treatmentRFC=1" = "solid"),
    labels = c("treatmentRFC=0" = "FC",     "treatmentRFC=1" = "RFC")
  ) +
  guides(
    color    = guide_legend(order = 1, override.aes = list(linetype = "solid")),
    linetype = guide_legend(order = 2, override.aes = list(color = "black"))
  )


# ==============================================================================
# SECTION: M-spline model with external aggregated data
# ==============================================================================

#' External data anchors the extrapolation by providing Binomial counts
#' r_j ~ Binomial(pi_j, n_j) over future time intervals.
#' Here the data are hypothetical (could come from registry, expert elicitation
#' or real-world evidence), representing survival in the FC arm beyond the trial.
#'
#' The model links these to the survival curve via:
#'   pi_j = S(t_j^stop | theta) / S(t_j^start | theta)
#' This is handled automatically inside survextrap when external= is supplied.

extdat = tibble(
  start     = c(100, 120, 150),
  stop      = c(120, 150, 180),
  n         = c(100, 100, 100),     # individuals at risk at the start of each interval
  r         = c(8,   6,   2),       # survivors at the end of each interval
  treatment = as.factor(c("FC", "FC", "FC"))
)

#' Re-run the M-spline with the external data; all other arguments unchanged.
msp_ext = survextrap(
  Surv(time, status) ~ treatment,
  data     = data,
  mspline  = mspline,
  chains   = 2,
  iter     = 4000,
  external = extdat
)

#' Extract posterior survival summaries from the external-data model.
S_ext = survextrap::survival(
  msp_ext,
  summ_fns = list(mean = mean, ~quantile(.x, probs = c(0.025, 0.975)))
)

#' Overlay the two M-spline variants (with and without external data).
#' The external-data model is pulled towards 0 at later times and shows
#' lower uncertainty in the extrapolation period.
p_base = plot(m, mods = c(1, 2), add.km = TRUE, t = seq(0, 180))

p2 = p_base +
  geom_line(
    data = S_ext |>
      mutate(
        model  = "M-spline (external evidence)",
        strata = case_when(
          treatment == "FC" ~ "treatmentRFC=0", TRUE ~ "treatmentRFC=1"
        ) |> factor(levels = c("treatmentRFC=1", "treatmentRFC=0"))
      ),
    aes(t, mean, col = model, linetype = strata), linewidth = 0.9
  ) +
  geom_line(
    data = S |>
      mutate(
        model  = "M-spline",
        strata = case_when(
          treatment == "FC" ~ "treatmentRFC=0", TRUE ~ "treatmentRFC=1"
        ) |> factor(levels = c("treatmentRFC=1", "treatmentRFC=0"))
      ),
    aes(t, mean, col = model, linetype = strata), linewidth = 0.9
  )

p2$layers[[1]] = NULL   # remove the parametric model layer for clarity
p2 +
  scale_linetype_manual(
    values = c("treatmentRFC=0" = "dashed", "treatmentRFC=1" = "solid"),
    labels = c("treatmentRFC=0" = "FC",     "treatmentRFC=1" = "RFC")
  ) +
  scale_color_manual(values = c("#1F77B4", "#FF7F0E"))

#' Posterior uncertainty for the external-data M-spline alone,
#' shown as a shaded ribbon.
S_ext |>
  mutate(
    model  = "M-spline (external evidence)",
    strata = case_when(
      treatment == "FC" ~ "treatmentRFC=0", TRUE ~ "treatmentRFC=1"
    ) |> factor(levels = c("treatmentRFC=1", "treatmentRFC=0"))
  ) |>
  ggplot(aes(t, mean, col = model, linetype = strata)) +
  geom_line(size = 0.9) +
  geom_ribbon(
    aes(ymin = `2.5%`, ymax = `97.5%`, group = strata),
    alpha = 0.2, linetype = 0, col = NA, show.legend = FALSE
  ) +
  theme_survHE() +
  theme(legend.position = c(0.55, 0.95), legend.justification = c("left", "top")) +
  labs(color = "Model", linetype = "Profile", x = "Time", y = "Survival") +
  scale_linetype_manual(
    values = c("treatmentRFC=0" = "dashed", "treatmentRFC=1" = "solid"),
    labels = c("treatmentRFC=0" = "FC",     "treatmentRFC=1" = "RFC")
  ) +
  guides(
    color    = guide_legend(order = 1, override.aes = list(linetype = "solid")),
    linetype = guide_legend(order = 2, override.aes = list(color = "black"))
  )
