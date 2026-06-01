#' ---
#' title: "Value of Information"
#' desc:  "Value of Information (VoI) analysis for HTA decision-making.
#'         Covers: EVPI computation from PSA simulations (Opportunity Loss);
#'         the Vaccine dataset EVPI curve via BCEA::evi.plot(); EVPPI via
#'         BCEA::evppi() and the Info Rank plot; a BART illustration for
#'         non-linear regression (dbarts); EVSI and ENBS using the voi
#'         package with the built-in Chemotherapy example; and standard
#'         power-based vs ENBS-based sample size calculations."
#' ---

library(tidyverse)
library(BCEA)
library(voi)


# ==============================================================================
# SECTION: PSA simulations and EVPI building blocks
# ==============================================================================

#' Simulates PSA output for a two-treatment decision model.
#' nb0 = net benefit for the status quo (t=1)
#' nb1 = net benefit for the new intervention (t=2)
#'
#' Population size (9685) converts per-person utility to aggregate MNB.
#' The random seed ensures the PSA table in the chapter is reproducible.

set.seed(140873)
nb0 = 9685 * round(rnorm(10000, 7.5, 3))   # NB for t=1 across 10,000 PSA sims
nb1 = 9685 * round(rnorm(10000, 8.0, 3))   # NB for t=2
ib  = nb1 - nb0                             # incremental benefit

#' Row-wise maximum NB (the "perfect information" decision for each simulation)
#' maxnb[s] = max(nb0[s], nb1[s]): the value a perfect-information decision
#' maker would achieve at simulation s.
maxnb = pmax(nb1, nb0)

#' Opportunity Loss (OL): how much we lose by not having perfect information.
#' OL[s] = maxnb[s] - nb_{t*}[s], where t* = arg max E[nb_t].
#' When t*=2 (the most cost-effective treatment on average), OL = maxnb - nb1.
#' OL is non-negative by construction and equals 0 whenever nb1 >= nb0.
ol = pmax(nb1, nb0) - nb1

#' EVPI = E[maxnb] - max(E[nb0], E[nb1]) = E[OL]
#' Three equivalent computations:
EVPI_from_maxnb = mean(maxnb) - max(mean(nb1), mean(nb0))
EVPI_from_OL    = mean(ol)
cat("EVPI (from maxnb):", format(EVPI_from_maxnb, digits = 2, nsmall = 2), "\n")
cat("EVPI (from OL):   ", format(EVPI_from_OL,    digits = 2, nsmall = 2), "\n")
cat("E[NB_1] =", format(mean(nb0), digits = 2), "\n")
cat("E[NB_2] =", format(mean(nb1), digits = 2), "\n")


# ==============================================================================
# SECTION: PSA summary table (VoI perspective)
# ==============================================================================

#' Constructs the PSA table adding the VoI columns: NB*(theta) and OL(theta).
#' Italics flag which arm has the higher NB row-by-row (to highlight OL).

show.col = 4
tab = tibble(
  iter  = c(format(seq(1, show.col)), NA, "S=1000"),
  pi0   = c(rbeta(show.col, 1, 2) |> round(2), NA, rbeta(1, 1, 2) |> round(2)),
  rho   = c(rbeta(show.col, 1, 2) |> round(2), NA, rbeta(1, 1, 2) |> round(2)),
  dots  = rep("...", show.col + 2),
  gamma = c(rbeta(show.col, 1, 2) |> round(2), NA, rbeta(1, 1, 2) |> round(2)),
  nb0   = c(nb0[1:show.col], NA, nb0[10000]),
  nb1   = c(nb1[1:show.col], NA, nb1[10000])
) |>
  mutate(
    maxnb = pmax(nb1, nb0),
    ol    = pmax(nb1, nb0) - nb1
  )

#' Flags for italic styling: row-wise indicator of which arm is optimal.
italics1 = tab$nb1 > tab$nb0; italics1[is.na(italics1)] = FALSE
italics0 = tab$nb0 > tab$nb1; italics0[is.na(italics0)] = FALSE

#' Adds an "Average" footer row with E[NB_0], E[NB_1], E[NB*] and E[OL].
tab2 = tab |> mutate(across(everything(), ~as.character(.)))
tab2 = tab2 |> tibble::add_row(
  iter  = "",
  pi0   = "", rho = "", dots = "Average:", gamma = "",
  nb0   = format(mean(nb0),   digits = 2, nsmall = 2),
  nb1   = format(mean(nb1),   digits = 2, nsmall = 2),
  maxnb = format(mean(maxnb), digits = 2, nsmall = 2),
  ol    = format(mean(ol),    digits = 2, nsmall = 2)
)

colnames(tab2) = c(
  "Sims", "theta_1", "theta_2", "...", "theta_Q",
  "NB_1(theta)", "NB_2(theta)", "NB*(theta)", "OL(theta)"
)

tab2 |>
  tinytable::tt() |>
  tinytable::group_tt(
    j = list("Parameter simulations" = 2:5,
             "Expected utility"       = 6:7,
             "VoI analysis"           = 8:9)
  ) |>
  tinytable::style_tt(i = which(italics0), j = 6, italic = TRUE) |>
  tinytable::style_tt(i = which(italics1), j = 7, italic = TRUE) |>
  tinytable::format_tt(replace = "...") |>
  tinytable::style_tt(i = 7, j = 4, colspan = 2, align = "r", bold = TRUE) |>
  tinytable::style_tt(j = 2:8, align = "c") |>
  tinytable::style_tt(i = 6, j = 1:9, line = "b")


# ==============================================================================
# SECTION: EVPI curve — Vaccine dataset (BCEA)
# ==============================================================================

#' Uses the Vaccine dataset bundled with BCEA to produce the EVPI curve.
#' bcea() creates the cost-effectiveness analysis object; evi.plot() plots
#' the EVPI as a function of the willingness-to-pay threshold k.
#'
#' The "kink" in the curve at k = m$kstar is the break-even point (BEP):
#' the value of k at which the optimal decision switches.  The two-segment
#' shape reflects that uncertainty costs more just before the BEP and less
#' after it, once one arm dominates more decisively.

data(Vaccine)
m = bcea(eff, cost, ref = 2)
evi.plot(m, graph = "gg")

#' The kink (BEP) and the EVPI scale.
cat("Break-even point k* =", m$kstar, "\n")
cat("Max EVPI ≈ £2.5 per person — indicating limited decision uncertainty.\n")


# ==============================================================================
# SECTION: EVPPI — regression-based approximation (BCEA + voi)
# ==============================================================================

#' Continues with the Vaccine model.
#' createInputs() extracts PSA simulations for all model parameters and
#' returns a list with:
#'   $mat        — S × Q matrix of parameter simulations
#'   $parameters — parameter names (for subsetting)
#'
#' evppi() fits a regression (GAM, GP or BART, controlled by the method=
#' argument) to approximate E[psi | phi] for the chosen subset of parameters
#' and returns the EVPPI as a function of k.

inp = createInputs(vaccine_mat, print_is_linear_comb = FALSE)

#' Compute the EVPPI for three focal parameters (beta.6, gamma.2, Trt.2.2.2).
#' The GAM/regression approximation is much cheaper than nested MC because it
#' only uses the existing S PSA simulations, not S_phi * S_psi new runs.
ev = BCEA::evppi(
  m,
  inp$parameters[grepl("beta.6.|gamma.2.|Trt.2.2.2.", inp$parameters)],
  inp$mat
)

#' Overlay EVPI (upper bound) and EVPPI on the same plot.
tibble(k = ev$k, voi = ev$evi,   type = "EVPI")  |>
  bind_rows(
    tibble(k = ev$k, voi = ev$evppi, type = "EVPPI")
  ) |>
  ggplot(aes(k, voi, linetype = type)) +
  geom_line() +
  xlab("Willingness to pay") + ylab("Expected VoI") +
  theme(legend.position = "bottom") + labs(linetype = "")

#' Info Rank plot: single-parameter EVPPI for all model parameters.
#' Bars ranked by EVPPI at the chosen k value — a guide to research priorities.
#' Note: VoI is non-additive, so the sum of single-parameter EVPPIs can exceed
#' the joint EVPPI for multiple parameters, or be lower if parameters interact.
info.rank(m, inp, graph = "gg")


# ==============================================================================
# SECTION: BART illustration — non-linear regression for EVPPI
# ==============================================================================

#' Illustrates how BART can recover a non-linear (sinusoidal) relationship.
#' In the EVPPI context, x ~ phi (focal parameters) and y ~ NB_t(theta);
#' the fitted BART ensemble approximates E[NB_t | phi], which is the
#' function g(phi) in eq-non-linear-evppi-regression.
#'
#' A single tree fits a piecewise constant function (panel b).
#' The ensemble of 200 trees (grey curve) recovers the smooth sinusoid.
#' BART uses regularisation priors to prevent any single tree from overfitting.

set.seed(1234)
c1 = 1.5; c2 = 3.5
x  = seq(.75, 5.75, .01)
y  = sin(x) + rnorm(length(x), mean = 0, sd = .2)

#' Fit BART with default settings (200 trees, 1000 posterior draws).
bart_model = dbarts::bart(
  x.train  = as.matrix(x),
  y.train  = y,
  keeptrees = TRUE
)
#' Posterior predictive mean across all trees and MCMC draws.
bart_preds = colMeans(predict(bart_model, as.matrix(x)))

#' Panel (a): binary decision tree diagram with two cutoffs c1, c2.
#' The tree assigns group means mu_1, mu_2, mu_3 to three x-regions.
ggplot() +
  annotate("label", x = 1, y = 2, label = "x < c1",
           size = 5.5, label.padding = unit(.35, "lines")) +
  annotate("label", x = 0, y = 1, label = "mu_1",
           size = 5.5, label.padding = unit(.35, "lines")) +
  annotate("label", x = 2, y = 1, label = "x < c2",
           size = 5.5, label.padding = unit(.35, "lines")) +
  annotate("label", x = 1, y = 0, label = "mu_2",
           size = 5.5, label.padding = unit(.35, "lines")) +
  annotate("label", x = 3, y = 0, label = "mu_3",
           size = 5.5, label.padding = unit(.35, "lines")) +
  xlim(-.5, 3.5) +
  geom_segment(aes(x = .85, y = 1.9,  xend = 0,  yend = 1.1)) +
  geom_segment(aes(x = 1.15, y = 1.9, xend = 2,  yend = 1.1)) +
  geom_segment(aes(x = 1.85, y = .9,  xend = 1,  yend = .1)) +
  geom_segment(aes(x = 2.15, y = .9,  xend = 3,  yend = .1)) +
  geom_label(aes(x = .25,  y = 1.5, label = "yes"), fill = "gray80", size = 4) +
  geom_label(aes(x = 1.75, y = 1.5, label = "no"),  fill = "gray80", size = 4) +
  geom_label(aes(x = 1.25, y = .5,  label = "yes"), fill = "gray80", size = 4) +
  geom_label(aes(x = 2.75, y = .5,  label = "no"),  fill = "gray80", size = 4) +
  theme_classic() +
  theme(
    axis.line.x = element_blank(), axis.line.y = element_blank(),
    axis.ticks  = element_blank(), axis.text   = element_blank()
  ) + xlab("") + ylab("")

#' Panel (b): observed data, single-tree piecewise fit, and BART ensemble.
tibble(x, y, pred = bart_preds) |>
  ggplot(aes(x, y)) +
  geom_point(col = "grey90") +
  geom_vline(xintercept = c1, linetype = 2) +
  geom_vline(xintercept = c2, linetype = 2) +
  #' Single-tree piecewise constant fit (group means)
  geom_segment(aes(x = -Inf, xend = c1,
                   y = mean(y[x < c1]),   yend = mean(y[x < c1])),
               linewidth = 1.4, col = "blue") +
  geom_segment(aes(x = c1, xend = c2,
                   y = mean(y[x > c1 & x < c2]), yend = mean(y[x > c1 & x < c2])),
               linewidth = 1.4, col = "blue") +
  geom_segment(aes(x = c2, xend = Inf,
                   y = mean(y[x > c2]),   yend = mean(y[x > c2])),
               linewidth = 1.4, col = "blue") +
  annotate("text", x = 1,           y = mean(y[x < c1]),
           label = "mu_1", vjust = -1.5, hjust = 0, size = 6, col = "blue") +
  annotate("text", x = (c2 + c1)/2, y = mean(y[x > c1 & x < c2]),
           label = "mu_2", vjust = -1.5, hjust = 0, size = 6, col = "blue") +
  annotate("text", x = (6 + c2)/2,  y = mean(y[x > c2]),
           label = "mu_3", vjust = -1.5, hjust = 0, size = 6, col = "blue") +
  #' BART ensemble prediction (smooth curve)
  geom_line(aes(x, pred), col = "grey45", linewidth = .75) +
  xlab("x") + ylab("y") +
  scale_x_continuous(breaks = c(c1, c2), labels = c("c1", "c2"))


# ==============================================================================
# SECTION: EVPPI simulation table (regression framing)
# ==============================================================================

#' Constructs the tbl-psa-evppi demonstration table:
#' the NB for a treatment plays the role of the regression "outcome" y,
#' while the PSA draws for the focal parameters phi_1, ..., phi_Q play the
#' role of regression "covariates" X.  The nuisance parameters psi_1, ...
#' are included to show the full data structure but are not used in the
#' regression.

show.col = 4
nbt = 6938 * round(rnorm(10000, 8, 3))

tab3 = tibble(
  iter  = c(format(seq(1, show.col)), NA, "S"),
  nbt   = c(nbt[1:show.col], NA, nbt[10000]),
  phi1  = c(rbeta(show.col, 1, 2)  |> round(2), NA, rbeta(1, 1, 2)  |> round(2)),
  phi2  = c(rgamma(show.col, 2, 5) |> round(2), NA, rgamma(1, 2, 5) |> round(2)),
  dots1 = rep("...", show.col + 2),
  phiQ  = c(rbeta(show.col, 1, 2)  |> round(2), NA, rbeta(1, 1, 2)  |> round(2)),
  psi1  = c(rnorm(show.col, 0, 2)  |> round(2), NA, rnorm(1, 0, 2)  |> round(2)),
  psi2  = c(rbeta(show.col, 1, 2)  |> round(2), NA, rbeta(1, 1, 2)  |> round(2)),
  dots2 = rep("...", show.col + 2),
  psiQ  = c(rbeta(show.col, 1, 2)  |> round(2), NA, rbeta(1, 1, 2)  |> round(2))
)

colnames(tab3) = c(
  "Sims", "NB_t(theta)",
  "phi_1", "phi_2", "...", "phi_Q_phi",
  "psi_1", "psi_2", "...", "psi_Q_psi"
)

tab3 |>
  tinytable::tt() |>
  tinytable::group_tt(
    j = list(
      "Outcome (y)"         = 2,
      "Covariates (X)"      = 3:6,
      "Nuisance parameters" = 7:10
    )
  ) |>
  tinytable::format_tt(replace = "...") |>
  tinytable::style_tt(j = 2:10, align = "c")


# ==============================================================================
# SECTION: EVSI — Chemotherapy example (voi package)
# ==============================================================================

#' Uses the Chemotherapy example from the voi package (vignette:
#' https://cran.r-project.org/web/packages/voi/vignettes/plots.html).
#'
#' chemo_evsi_or: pre-computed EVSI object for a two-arm trial with binary
#'   outcome (side-effects), varying sample size n and WTP k.
#'   The proposed study updates only the log OR of side effects between arms.
#' chemo_cea_501: cost-effectiveness outputs with 501 PSA samples.
#'
#' The EVPI acts as the upper bound for the EVPPI, which in turn bounds
#' the EVSI from above (for a trial of finite size).

evpi_df  = evpi(outputs = chemo_cea_501)
evppi_df = attr(chemo_evsi_or, "evppi")

#' Panel 1: EVSI by WTP threshold for different sample sizes.
#' The EVSI peaks near the break-even point (~k=20000) and flattens for
#' larger k where one treatment dominates more decisively.
ggplot(chemo_evsi_or, aes(x = k, y = evsi, group = n, color = n)) +
  geom_line() +
  scale_colour_gradient(low = "skyblue", high = "darkblue") +
  xlab("Willingness to pay") + ylab("EVSI per person") +
  geom_line(data = evpi_df,  aes(x = k, y = evpi),  color = "black",   linewidth = 1.2, inherit.aes = FALSE) +
  geom_line(data = evppi_df, aes(x = k, y = evppi), color = "darkblue", linewidth = 1.2, inherit.aes = FALSE) +
  labs(color = "Sample size") + xlim(0, 54000) +
  annotate("text", x = 50000, y = 125, label = "EVPI",  hjust = 0) +
  annotate("text", x = 50000, y = 107, label = "EVPPI", color = "darkblue", hjust = 0) +
  annotate("text", x = 50000, y =  50, label = "EVSI",  color = "darkblue", hjust = 0) +
  theme(
    legend.position.inside = c(.85, .8), legend.position = "inside",
    legend.background      = element_blank(),
    legend.text  = element_text(size = 7),
    legend.title = element_text(size = 8)
  )

#' Panel 2: EVSI by sample size for different WTP thresholds.
#' The curves flatten around n=1500: diminishing marginal returns from
#' larger samples, since the parameter of interest becomes well-estimated.
chemo_evsi_or |>
  dplyr::filter(k %in% seq(15000, 50000, by = 5000)) |>
  ggplot(aes(x = n, y = evsi, group = k, color = k)) +
  geom_line() +
  scale_colour_gradient(low = "skyblue", high = "darkblue",
                        breaks = c(15000, 50000)) +
  xlab("Sample size") + ylab("EVSI per person") +
  labs(color = "Willingness-to-pay") + ylim(0, 450) +
  theme(
    legend.position.inside = c(.12, .8), legend.position = "inside",
    legend.background      = element_blank(),
    legend.text  = element_text(size = 7),
    legend.title = element_text(size = 8)
  )


# ==============================================================================
# SECTION: Standard power-based sample size vs ENBS-based sample size
# ==============================================================================

#' Standard frequentist power calculation for a two-arm trial with binary
#' outcome.  Assumes a background side-effect rate of theta1=0.47 and targets
#' an OR=0.55 reduction under chemotherapy.

theta1 = 0.47
OR     = 0.55
theta2 = (theta1 * OR) / (1 - theta1 + (theta1 * OR))

#' power.prop.test() computes power for a two-proportion z-test.
#' n=seq(0,1500,10) returns the power curve over a grid of sample sizes.
ssc = power.prop.test(
  p1        = theta1,
  p2        = theta2,
  sig.level = 0.05,
  n         = seq(0, 1500, 10)
)

#' Power curve: the 0.8 threshold corresponds to ~190 individuals per arm.
tibble(n = ssc$n, power = ssc$power) |>
  ggplot(aes(n, power)) +
  geom_line() +
  xlab("Proposed sample size") + ylab("Power") +
  geom_hline(yintercept = .8, linetype = 2) +
  annotate(
    "segment",
    x    = ssc$n[which(ssc$power >= .8) |> min()],
    xend = ssc$n[which(ssc$power >= .8) |> min()],
    y    = -Inf, yend = 0.8, linetype = "dashed"
  ) +
  scale_x_continuous(
    breaks = c(pretty(ssc$n), ssc$n[which(ssc$power >= .8) |> min()]),
    labels = c(pretty(ssc$n), ssc$n[which(ssc$power >= .8) |> min()])
  ) +
  scale_y_continuous(
    breaks = c(pretty(ssc$power), 0.8),
    labels = c(pretty(ssc$power), 0.8)
  )

cat("Power-based n per arm:", ssc$n[which(ssc$power >= .8) |> min()], "\n")

# ==============================================================================
# SECTION: Expected Net Benefit of Sampling (ENBS)
# ==============================================================================

#' ENBS = EVSI - expected cost of the study.
#' enbs() takes the EVSI object, adds study cost assumptions and scales to
#' the target population and time horizon.
#'
#' Cost assumptions (from the voi vignette):
#'   Fixed setup cost:   £5M – £10M
#'   Cost per patient:   £28K – £42K  (the chapter states £14K – £21K per arm;
#'                                      total per patient = 2 arms × these values)
#'   Target population:  46,000 patients
#'   Time horizon:       10 years (3.5% annual discount rate built in)

nbs = enbs(
  chemo_evsi_or,
  costs_setup = c(5000000, 10000000),
  costs_pp    = c(28000, 42000),
  pop         = 46000,
  time        = 10
)

#' Extract the ENBS at k=20000 and compute credible intervals from the
#' Normal approximation (mean = enbs, sd = sd of the ENBS estimate).
nbs |>
  dplyr::filter(k == 20000) |>
  mutate(
    q975 = qnorm(0.975, enbs, sd),
    q75  = qnorm(0.75,  enbs, sd),
    q25  = qnorm(0.25,  enbs, sd),
    q025 = qnorm(0.025, enbs, sd)
  ) |>
  ggplot(aes(y = enbs, x = n)) +
  geom_ribbon(aes(ymin = q025, ymax = q975), fill = "gray80", alpha = .45) +
  geom_ribbon(aes(ymin = q25,  ymax = q75),  fill = "gray50", alpha = .45) +
  geom_line() +
  xlab("Sample size") + ylab("Expected net benefit of sampling") +
  scale_y_continuous(labels = scales::dollar_format(prefix = "£")) +
  annotate("text", x = 1100, y = 85000000, label = "95% credible interval") +
  annotate("text", x = 1250, y = 70000000, label = "50% credible interval")

#' ENBS at the power-based sample size (n ≈ 190 per arm).
enbsssc = nbs |> dplyr::filter(k == 20000, n == 200) |> pull(enbs)

#' Optimal n and maximum ENBS.
nopt = (nbs |> dplyr::filter(k == 20000))$n[
  which((nbs |> dplyr::filter(k == 20000) |> as_tibble() |> pull(enbs)) ==
        (nbs |> dplyr::filter(k == 20000) |> as_tibble() |> pull(enbs) |> max()))
]
enbsopt = (nbs |> dplyr::filter(k == 20000))$enbs[
  which((nbs |> dplyr::filter(k == 20000) |> as_tibble() |> pull(enbs)) ==
        (nbs |> dplyr::filter(k == 20000) |> as_tibble() |> pull(enbs) |> max()))
]

cat("ENBS at power-based n=200 per arm: £", format(enbsssc, big.mark=","), "\n")
cat("Optimal n:", nopt, "per arm\n")
cat("Maximum ENBS: £", format(enbsopt, digits=1, nsmall=0, big.mark=","), "\n")
