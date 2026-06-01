#' ---
#' title: "Population adjustment"
#' desc:  "Methods for indirect treatment comparisons under effect-modifier
#'         imbalance across trial populations. Code covers: (1) the
#'         collapsibility example illustrating the distinction between
#'         conditional and marginal treatment effects on the RD, RR and OR
#'         scales; (2) Monte Carlo vs Quasi-Monte Carlo integration
#'         illustrated on a 2D indicator function (Halton sequence); and
#'         (3) derivation of the population-average conditional and marginal
#'         estimands from the collapsibility example, using bmhe::logit()."
#' ---

library(tidyverse)
rm(list=ls()) # resets workspace to remove existing variables


# ==============================================================================
# SECTION: Collapsibility — contingency table setup
# ==============================================================================

#' Fictional two-arm trial (Treatment vs Placebo) with smoking as a covariate.
#' The table records the number of successes and failures in each
#' treatment × smoking stratum.
#'
#' Key findings illustrated:
#'   - Smoking is a *prognostic factor* (affects absolute outcomes) but under
#'     the data as given, treatment is equally distributed across strata, so
#'     smoking is *not* a confounder.
#'   - On the Risk Difference (RD) and Odds Ratio (OR) scales, smoking is
#'     *not* an effect modifier (stratum-specific effects are equal).
#'   - On the Risk Ratio (RR) scale, smoking IS an effect modifier (stratum-
#'     specific RRs differ).
#'   - The RD is *collapsible*: the crude overall RD equals the weighted
#'     average of the stratum-specific RDs.
#'   - The OR is *not collapsible*: the crude OR differs from the weighted
#'     average, even in the absence of confounding or effect modification.

tab = tibble(
  trt               = c("Treatment", "Placebo"),
  rsp_non_smokers   = c(135,  50),
  nrsp_non_smokers  = c( 15,  50),
  rsp_smokers       = c( 45,   6),
  nrsp_smokers      = c( 45,  54)
) |>
  mutate(
    rsp  = rsp_non_smokers  + rsp_smokers,
    nrsp = nrsp_non_smokers + nrsp_smokers
  )

#' Display table with row totals and grouped column headers.
tab |>
  add_row(
    trt               = "Total",
    rsp_non_smokers   = sum(tab$rsp_non_smokers),
    nrsp_non_smokers  = sum(tab$nrsp_non_smokers),
    rsp_smokers       = sum(tab$rsp_smokers),
    nrsp_smokers      = sum(tab$nrsp_smokers),
    rsp               = sum(tab$rsp),
    nrsp              = sum(tab$nrsp)
  ) |>
  tinytable::tt() |>
  tinytable::group_tt(j = list("Non smokers" = 2:3, "Smokers" = 4:5, "Total" = 6:7)) |>
  setNames(c("", rep(c("Success", "Failure"), 3))) |>
  tinytable::style_tt(i = 3, line = "t")


# ==============================================================================
# SECTION: Collapsibility — stratum-specific effect measures
# ==============================================================================

#' Column index shortcuts (for clarity in the calculations below):
#'   tab[, 2] = rsp_non_smokers
#'   tab[, 3] = nrsp_non_smokers
#'   tab[, 4] = rsp_smokers
#'   tab[, 5] = nrsp_smokers
#'   tab[, 6] = rsp (total)
#'   tab[, 7] = nrsp (total)
#'
#' Row 1 = Treatment arm; Row 2 = Placebo arm.

# --- Stratum sizes ---
n_ns    = sum(tab[, 2:3])   # total in non-smoker stratum
n_s     = sum(tab[, 4:5])   # total in smoker stratum
n_total = sum(tab[, 2:5])   # overall total

# --- Stratum weights ---
w_ns = n_ns / n_total
w_s  = n_s  / n_total

# --- Success rates per stratum and arm ---
p_trt_ns = tab[1, 2] / sum(tab[1, 2:3])   # Treatment, non-smokers
p_pbo_ns = tab[2, 2] / sum(tab[2, 2:3])   # Placebo, non-smokers
p_trt_s  = tab[1, 4] / sum(tab[1, 4:5])   # Treatment, smokers
p_pbo_s  = tab[2, 4] / sum(tab[2, 4:5])   # Placebo, smokers

# --- Risk Difference (RD) ---
#' RD is collapsible: crude = weighted average of strata-specific.
RD_ns    = (p_trt_ns - p_pbo_ns) |> round(1)
RD_s     = (p_trt_s  - p_pbo_s)  |> round(1)
RD_crude = (tab[1, 6] / sum(tab[1, 6:7]) - tab[2, 6] / sum(tab[2, 6:7])) |> round(1)
RD_wtd   = (w_ns * RD_ns + w_s * RD_s) |> round(1)

cat(paste("RD: stratum NS =", RD_ns, ", stratum S =", RD_s, "\n"))
cat(paste("RD crude =", RD_crude, " | weighted average =", RD_wtd, "\n"))
cat(paste("Collapsible:", isTRUE(all.equal(as.numeric(RD_crude), as.numeric(RD_wtd))), "\n"))

# --- Odds Ratio (OR) ---
#' OR is NOT collapsible: crude differs from weighted average on log-OR scale,
#' even though smoking is not an effect modifier on the OR scale.
OR_ns    = (tab[1, 2] * tab[2, 3]) / (tab[1, 3] * tab[2, 2])
OR_s     = (tab[1, 4] * tab[2, 5]) / (tab[1, 5] * tab[2, 4])
OR_crude = (tab[1, 6] * tab[2, 7]) / (tab[1, 7] * tab[2, 6]) |> round(1)
#' Weighted average on log scale, then exponentiated back
OR_wtd   = exp(w_ns * log(OR_ns) + w_s * log(OR_s))

cat(paste("\nOR: stratum NS =", OR_ns, ", stratum S =", OR_s, "\n"))
cat(paste("OR crude =", OR_crude, " | weighted average =", round(OR_wtd, 4), "\n"))
cat(paste("Collapsible:", isTRUE(all.equal(as.numeric(OR_crude), round(OR_wtd, 1))), "\n"))

# --- Risk Ratio (RR) ---
#' RR: smoking IS an effect modifier on this scale (stratum-specific RRs differ).
RR_ns = (p_trt_ns / p_pbo_ns) |> round(1)
RR_s  = (p_trt_s  / p_pbo_s)  |> round(1)
cat(paste("\nRR: stratum NS =", RR_ns, ", stratum S =", RR_s,
    " — effect modification present on this scale\n"))


# ==============================================================================
# SECTION: Population-average conditional and marginal estimands (ML-NMR)
# ==============================================================================

#' Continuing from the collapsibility example, computes the two key estimands
#' that ML-NMR and g-computation target, respectively:
#'
#' (a) Population-average CONDITIONAL treatment effect (d_BA):
#'     d_BA = sum_x [logit(mu_B(x)) - logit(mu_A(x))] * f(x)
#'     This is the average of the individual-level logit contrasts, weighted
#'     by the target covariate distribution.  With a logit link, this is
#'     collapsible iff the response surface is linear in x on the logit scale.
#'
#' (b) Population-average MARGINAL treatment effect (delta_BA):
#'     first compute marginal means: mu_bar_t = sum_x mu_t(x) * f(x)
#'     then apply the link: delta_BA = g(mu_bar_B) - g(mu_bar_A)
#'     This is the contrast of the overall (pooled) marginal means on the
#'     link scale and is always the correct target for economic models.
#'
#' Because the OR/logit is non-collapsible, d_BA != delta_BA.

mu_B_NS = (tab[1, 2] / sum(tab[1, 2:3])) |> round(1)
mu_A_NS = (tab[2, 2] / sum(tab[2, 2:3])) |> round(1)
mu_B_S  = (tab[1, 4] / sum(tab[1, 4:5])) |> round(1)
mu_A_S  = (tab[2, 4] / sum(tab[2, 4:5])) |> round(1)

# Covariate distribution in the target population
f_ns = (sum(tab[, 2:3]) / sum(tab[, 2:5])) |> round(3)
f_s  = (sum(tab[, 4:5]) / sum(tab[, 2:5])) |> round(3)

#' (a) Conditional estimand — weighted sum of logit contrasts
d_BA = (
  (bmhe::logit(mu_B_NS) - bmhe::logit(mu_A_NS)) * f_ns +
  (bmhe::logit(mu_B_S)  - bmhe::logit(mu_A_S))  * f_s
) |> round(4)
cat(paste("Population-average conditional treatment effect (d_BA):", d_BA, "\n"))

#' (b) Marginal estimand — contrast of marginal means on link scale
bar_mu_B = (mu_B_NS * f_ns + mu_B_S * f_s) |> round(2)
bar_mu_A = (mu_A_NS * f_ns + mu_A_S * f_s) |> round(2)

delta_BA = (bmhe::logit(bar_mu_B) - bmhe::logit(bar_mu_A)) |> round(4)
cat(paste("Population-average marginal treatment effect (delta_BA):", delta_BA, "\n"))

#' Illustrates non-collapsibility: d_BA != delta_BA
cat(paste("Non-collapsibility: d_BA (", d_BA, ") != delta_BA (", delta_BA, ")\n"))

#' Marginal means also coincide with the crude "Total" column in the table.
cat(paste("bar_mu_B =", bar_mu_B,
    "(cf. crude:", round(tab[1, 6] / sum(tab[1, 6:7]), 2), ")\n"))
cat(paste("bar_mu_A =", bar_mu_A,
    "(cf. crude:", round(tab[2, 6] / sum(tab[2, 6:7]), 2), ")\n"))


# ==============================================================================
# SECTION: Monte Carlo vs Quasi-Monte Carlo integration
# ==============================================================================

#' Illustrates the efficiency advantage of QMC over MC for estimating a
#' two-dimensional integral using n=100 function evaluations.
#'
#' Target: integral of f(x,y) = 1{(x-0.5)^2 + (y-0.5)^2 <= 0.09} over [0,1]^2.
#' Analytic value: pi * 0.3^2 ≈ 0.2827.
#'
#' MC: points sampled uniformly at random (pseudo-random).
#' QMC: points from a Halton sequence (base 2 for x, base 3 for y).
#'   The Halton sequence is a low-discrepancy sequence that fills the unit
#'   square more evenly than pseudo-random numbers, which tend to cluster.
#'   Convergence rate: O(n^{-1}) for smooth integrands vs O(n^{-1/2}) for MC.

set.seed(413)

#' Indicator function: 1 inside circle of radius 0.3 centred at (0.5, 0.5).
f = function(x, y) {
  ifelse((x - 0.5)^2 + (y - 0.5)^2 <= 0.3^2, 1, 0)
}

n        = 100
true_val = pi * 0.3^2   # analytic value of the integral

# --- Monte Carlo: pseudo-random uniform points ---
mc_x        = runif(n)
mc_y        = runif(n)
mc_vals     = f(mc_x, mc_y)
mc_estimate = mean(mc_vals)

# --- Quasi-Monte Carlo: Halton sequence ---
#' halton(n, base) generates the first n terms of the Halton sequence in
#' the given base, exploiting the van der Corput construction: each integer
#' i is expressed in the given base and then reflected about the decimal point.
halton = function(n, base) {
  result = numeric(n)
  for (i in 1:n) {
    f_val = 1; r = 0; j = i
    while (j > 0) {
      f_val = f_val / base
      r     = r + f_val * (j %% base)
      j     = floor(j / base)
    }
    result[i] = r
  }
  result
}

qmc_x        = halton(n, 2)   # base 2 for x dimension
qmc_y        = halton(n, 3)   # base 3 for y dimension (coprime bases important)
qmc_vals     = f(qmc_x, qmc_y)
qmc_estimate = mean(qmc_vals)

cat("True value:", round(true_val, 4), "\n")
cat("MC estimate:", round(mc_estimate, 4),
    " | Error:", round(abs(mc_estimate - true_val), 4), "\n")
cat("QMC estimate:", round(qmc_estimate, 4),
    " | Error:", round(abs(qmc_estimate - true_val), 4), "\n")

# --- Circle boundary for plot overlay ---
theta_seq   = seq(0, 2 * pi, length.out = 100)
circle_data = data.frame(
  x = 0.5 + 0.3 * cos(theta_seq),
  y = 0.5 + 0.3 * sin(theta_seq)
)

# --- Data frames for plotting ---
mc_data = data.frame(
  x = mc_x, y = mc_y, value = mc_vals,
  point_type = ifelse(mc_vals == 1, "f(x,y) = 1", "f(x,y) = 0")
)
qmc_data = data.frame(
  x = qmc_x, y = qmc_y, value = qmc_vals,
  point_type = ifelse(qmc_vals == 1, "f(x,y) = 1", "f(x,y) = 0")
)

#' Panel 1: Monte Carlo — random points show visible clustering and gaps.
p1 = ggplot() +
  geom_polygon(data = circle_data, aes(x = x, y = y),
               fill = "lightblue", alpha = 0.3, color = "steelblue", linewidth = 0.5) +
  geom_point(data = mc_data, aes(x = x, y = y, shape = point_type),
             size = 1.5, alpha = 0.8) +
  scale_shape_discrete(name = NULL) +
  coord_fixed() +
  theme(
    legend.position = "bottom",
    legend.text  = element_text(size = 15),
    axis.title   = element_text(size = 14),
    axis.text    = element_text(size = 12)
  ) +
  labs(
    x     = "x", y = "y",
    title = sprintf(
      "Analytic (true) value: %.4f; Estimate: %.4f; Error: %.4f",
      true_val, mc_estimate, abs(mc_estimate - true_val)
    )
  )
p1

#' Panel 2: Quasi-Monte Carlo — Halton points cover space more evenly.
p2 = ggplot() +
  geom_polygon(data = circle_data, aes(x = x, y = y),
               fill = "lightblue", alpha = 0.3, color = "steelblue", linewidth = 0.5) +
  geom_point(data = qmc_data, aes(x = x, y = y, shape = point_type),
             size = 1.5, alpha = 0.8) +
  scale_shape_discrete(name = NULL) +
  coord_fixed() +
  theme(
    legend.position = "bottom",
    legend.text  = element_text(size = 15),
    axis.title   = element_text(size = 14),
    axis.text    = element_text(size = 12)
  ) +
  labs(
    x     = "x", y = "y",
    title = sprintf(
      "Analytic (true) value: %.4f; Estimate: %.4f; Error: %.4f",
      true_val, qmc_estimate, abs(qmc_estimate - true_val)
    )
  )
p2
