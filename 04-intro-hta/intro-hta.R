#' ---
#' title: "Introduction to Health Technology Assessment"
#' desc:  "Introduces the statistical and decision-theoretic framework for
#'         health technology assessment (HTA). Code covers: QALY computation
#'         from utility trajectories; PSA simulation scaffolding; cost-
#'         effectiveness plane and ICER visualisation; and CEAC, CEAF and CEEF
#'         plots using the BCEA package."
#' ---

library(tidyverse)
library(BCEA)


# ==============================================================================
# SECTION: QALY computation
# ==============================================================================

#' Computes Quality-Adjusted Life Years (QALYs) from a sequence of utility
#' measurements taken at irregular time points.
#'
#' QALYs are the standard measure of health benefit in HTA.  They capture both
#' the quantity and quality of life by treating each observed utility score as
#' the height of a trapezoid and the gap to the next measurement as its width.
#' The total QALYs equal the area under the piecewise-linear utility curve:
#'
#'   QALYs = sum_j  [(u_j + u_{j-1}) / 2]  *  delta_j
#'
#' where delta_j = time_{j} - time_{j-1} is the length of the j-th interval.

# Measurement time points (years) and corresponding utility scores (0-1 scale)
years = c(0, 1.2, 2.3, 4, 7.1, 7.1)
qols  = c(0.88, 0.71, 0.38, 0.22, 0, 0)

#' Plots the utility trajectory as a shaded polygon (the "area under the curve").
#' Braces annotate the interval width (delta_j) and the trapezoidal height
#' ((u_j + u_{j-1}) / 2) for one representative interval.
tibble(years = years, qols = qols) |>
  tibble::add_row(years = c(7.1, 0), qols = c(0, 0)) |>
  ggplot(aes(x = years, y = qols)) +
  geom_polygon(fill = "#8080B3CC") +
  geom_point(data = tibble(x = years[1:5]), y = qols[1:5], aes(x, y)) +
  geom_segment(aes(x = 0,     xend = years, y = qols,  yend = qols),  lty = 2, size = .5, col = "grey40") +
  geom_segment(aes(x = years, xend = years, y = 0,     yend = qols),  lty = 2, size = .5, col = "grey40") +
  xlim(0, 10) + ylim(0, 1) +
  xlab("Time (years)") + ylab("Quality of life (0-1 scale)") +
  ggbrace::stat_brace(
    data = tibble(x = c(4, 5.55, 7.1), y = c(0.22, .28, .22)),
    aes(x = x, y = y), outside = FALSE
  ) +
  ggbrace::stat_brace(
    data = tibble(x = c(2.3, 2.7, 2.3), y = c(.38, .4625, .545)),
    aes(x, y), outside = FALSE, rotate = 90
  ) +
  annotate("text", x = 5.55, y = .35,    label = "${\\delta_j}$") +
  annotate("text", x = 3.7,  y = .4625,  label = "$\\displaystyle{\\frac{u_{j}+u_{j-1}}{2}}$") +
  geom_segment(aes(x = 7.1, xend = 7.1, y = 0,   yend = .22), lty = 2, size = .5, col = "grey40") +
  geom_segment(aes(x = 0,   xend = 2.3, y = .545, yend = .545), lty = 2, size = .5, col = "grey40")

#' Computes the total QALYs numerically using the trapezoid rule.
#' lag() aligns each measurement with the previous one so both components of
#' each trapezoid can be computed in a single mutate() call.
qalys = tibble(years = years, qols = qols) |>
  mutate(
    time     = years - dplyr::lag(years, default = 0),
    qaly_ind = (qols + dplyr::lag(qols, default = 0)) / 2,
    qaly     = qaly_ind * time
  ) |>
  pull(qaly) |>
  sum()

qalys   # Total QALYs: the individual lives max(years) years but this
        # is equivalent to 'qalys' years in perfect health


# ==============================================================================
# SECTION: ICER and cost-effectiveness plane
# ==============================================================================

#' Draws the cost-effectiveness plane for four pairwise comparisons of a
#' reference intervention A against B, C, D and E.
#' The NE and SW quadrants are shaded grey to highlight the region where the
#' ICER does not determine dominance unambiguously.
tibble(
  x   = c(100,  200,  200, -100),
  y   = c(-100, -100, -200,  150),
  col = c("magenta", "red", "blue", "green4")
) |>
  ggplot(aes(x, y)) +
  scale_colour_identity() +
  theme(
    axis.line = element_blank(), axis.text.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks = element_blank(),
    legend.position = "none", panel.background = element_blank(),
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), plot.background = element_blank()
  ) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_segment(aes(x = 0,    y = -Inf, xend = 0,   yend = Inf),
               arrow = arrow(length = unit(0.30, "cm"), type = "closed"), col = "black") +
  geom_segment(aes(x = -Inf, y = 0,    xend = Inf, yend = 0),
               arrow = arrow(length = unit(0.30, "cm"), type = "closed"), col = "black") +
  geom_segment(aes(x = x, xend = x, y = 0, yend = y), linetype = "dashed", linewidth = .1, col = "black") +
  geom_segment(aes(x = 0, xend = x, y = y, yend = y), linetype = "dashed", linewidth = .1, col = "black") +
  geom_point(size = 3) +
  xlab("Effectiveness differential") + ylab("Cost differential") +
  annotate("text", -100,  150, label = "ICER[AE] *' = -1.5'", vjust = -0.75, parse = TRUE) +
  annotate("text",  100, -100, label = "ICER[AB] *' = -1'",   vjust =  1.75, parse = TRUE) +
  annotate("text",  200, -100, label = "ICER[AC] *' = -1/2'", hjust = -0.25, parse = TRUE) +
  annotate("text",  200, -200, label = "ICER[AD] *' = -1'",   hjust = -0.25, parse = TRUE) +
  annotate("text",  100, 0, label = "100",  vjust = -.8, size = 3) +
  annotate("text",  200, 0, label = "200",  vjust = -.8, size = 3) +
  annotate("text", -100, 0, label = "-100", vjust = 1.3, size = 3) +
  annotate("text", 0,  150, label = "150",  hjust = -.8, size = 3) +
  annotate("text", 0, -100, label = "-100", hjust = 1.3, size = 3) +
  annotate("text", 0, -200, label = "-200", hjust = 1.3, size = 3) +
  annotate("text",  120,  100, label = "+ Effectiveness but + Costs\n??") +
  annotate("text", -120, -100, label = "- Effectiveness but - Costs\n??") +
  geom_rect(aes(xmin = 0,    xmax = Inf,  ymin = 0,    ymax = Inf),  fill = "lightgrey", alpha = .05, color = NA) +
  geom_rect(aes(xmin = -Inf, xmax = 0,    ymin = -Inf, ymax = 0),    fill = "lightgrey", alpha = .05, color = NA) +
  xlim(-170, 290)


# ==============================================================================
# SECTION: PSA simulation scaffolding (illustrative)
# ==============================================================================

#' Simulates 10,000 PSA iterations for two interventions (t=1 and t=2).
#' Each iteration draws a value for each model parameter from its assumed
#' distribution and computes the expected net benefit for both interventions.
#'
#' In a real analysis the net benefit distributions would come directly from
#' MCMC posterior simulations; here they are stand-ins for illustration.
#' The population size (9685) converts per-person values to aggregate scale.
set.seed(140873)
nb0 = 9685 * round(rnorm(10000, 7.5, 3))  # NB for t=1 across 10,000 PSA sims
nb1 = 9685 * round(rnorm(10000, 8.0, 3))  # NB for t=2
ib  = nb1 - nb0                            # incremental benefit

#' Summarises the two expected net benefits and the EIB.
#' The intervention with the higher average NB is the cost-effective choice.
mean(nb0)  # E[NB_1]: expected net benefit for status quo
mean(nb1)  # E[NB_2]: expected net benefit for new intervention
mean(ib)   # EIB: positive => t=2 is cost-effective given current evidence

#' Constructs a display table of PSA iterations (tbl-psa in the chapter).
#' Shows a selection of rows alongside computed NB columns and averages.
show.col = 4
tab = tibble(
  iter  = c(format(seq(1, show.col)), NA, format(1000)),
  pi0   = c(rbeta(show.col, 1, 2) |> round(2), NA, rbeta(1, 1, 2) |> round(2)),
  rho   = c(rbeta(show.col, 1, 2) |> round(2), NA, rbeta(1, 1, 2) |> round(2)),
  dots  = rep("\\(\\ldots\\)", show.col + 2),
  gamma = c(rbeta(show.col, 1, 2) |> round(2), NA, rbeta(1, 1, 2) |> round(2)),
  nb0   = c(nb0[1:show.col], NA, nb0[10000]),
  nb1   = c(nb1[1:show.col], NA, nb1[10000])
) |>
  mutate(ib = c(ib[1:show.col], NA, ib[10000]))

#' Flags which intervention is optimal row-by-row (for italics in the table).
italics1 = tab$nb1 > tab$nb0; italics1[is.na(italics1)] = FALSE
italics0 = tab$nb0 > tab$nb1; italics0[is.na(italics0)] = FALSE

#' Appends a summary row with the overall expected net benefits and EIB.
tab2 = tab |> mutate(across(everything(), ~as.character(.)))
tab2 = tab2 |> add_row(
  iter  = "",
  pi0   = "", rho = "", dots = "Average:", gamma = "",
  nb0   = paste0("\\(\\mathcal{NB}_1\\!=\\,\\)", format(mean(nb0), digits = 2, nsmall = 2, big.mark = "")),
  nb1   = paste0("\\(\\mathcal{NB}_2\\!=\\,\\)", format(mean(nb1), digits = 2, nsmall = 2, big.mark = "")),
  ib    = paste0("\\(\\text{EIB}\\!=\\,\\)",     format(mean(ib),  digits = 2, nsmall = 2, big.mark = ""))
)

colnames(tab2) = c(
  "Sims", "\\(\\theta_1\\)", "\\(\\theta_2\\)", "\\(\\ldots\\)", "\\(\\theta_Q\\)",
  "\\(\\text{NB}_1(\\boldsymbol{\\theta})\\)",
  "\\(\\text{NB}_2(\\boldsymbol{\\theta})\\)",
  "\\(\\text{IB}(\\boldsymbol{\\theta})\\)"
)

#' Renders the table with tinytable, applying grouped column headers and
#' italic/bold styling to highlight which intervention is row-wise optimal.
tab2 |> tinytable::tt() |>
  tinytable::group_tt(j = list("Parameter simulations" = 2:5, "Expected utility" = 6:7)) |>
  tinytable::style_tt(j = 2:8, align = "c") |>
  tinytable::style_tt(i = which(italics0), j = 6, italic = TRUE) |>
  tinytable::style_tt(i = 7, j = 7, bold = TRUE) |>
  tinytable::style_tt(i = which(italics1), j = 7, italic = TRUE) |>
  tinytable::format_tt(replace = "\\(\\ldots\\)") |>
  tinytable::style_tt(i = 7, j = 4, colspan = 2, align = "r", bold = TRUE)


# ==============================================================================
# SECTION: Cost-effectiveness plane using BCEA (Vaccine dataset)
# ==============================================================================

#' Uses the built-in Vaccine dataset from BCEA to illustrate the
#' cost-effectiveness plane.
#'
#' bcea() takes matrices of effectiveness and cost simulations (rows = PSA
#' iterations, columns = interventions) and ref= specifies the reference
#' intervention index.  It returns a structured object that BCEA's plot
#' functions consume directly.

data(Vaccine)
m = bcea(eff = eff, cost = cost, ref = 2)

#' Extracts the incremental effectiveness and cost from the bcea object for
#' manual plotting of the CE plane.
df = data.frame(de = m$delta_e, dc = m$delta_c) |>
  as_tibble() |>
  rename(de = 1, dc = 2)

#' Panel 1: base-case ICER as a single point (point estimate only, no cloud).
p = ggplot(data = df, aes(de, dc)) +
  geom_point(col = "gray96", size = .4) +
  theme_classic() +
  theme(
    axis.line = element_blank(), axis.text.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks = element_blank(),
    legend.position = "none", panel.background = element_blank(),
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), plot.background = element_blank()
  ) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_segment(aes(x = 0,    y = -Inf, xend = 0,   yend = Inf),  arrow = arrow(length = unit(0.30, "cm"), type = "closed")) +
  geom_segment(aes(x = -Inf, y = 0,    xend = Inf, yend = 0),    arrow = arrow(length = unit(0.30, "cm"), type = "closed")) +
  geom_point(aes(x = mean(de), y = mean(dc)), col = "red", size = 2) +
  labs(title = "", y = "Cost differential", x = "Effectiveness differential") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15)
  ) +
  annotate("text", mean(df$de), mean(df$dc),
           label = "$\\mbox{ICER}=\\frac{ \\mbox{E}[\\Delta_c]}{ \\mbox{E}[\\Delta_e]}$",
           hjust = -.05, size = 6) +
  annotate("text", mean(df$de), mean(df$dc),
           label = "$\\phantom{\\mbox{ICER}}=\\mbox{Cost per outcome}$",
           hjust = -.05, vjust = 2.5, size = 6) +
  annotate("text", 0, Inf,  label = "$\\Delta_c$", hjust = 2,  vjust = 1, size = 5) +
  annotate("text", Inf, 0,  label = "$\\Delta_e$", hjust = 1,  vjust = 2, size = 5) +
  geom_segment(aes(x = mean(de), xend = mean(de), y = 0,         yend = mean(dc)), linetype = 2, size = .25) +
  geom_segment(aes(x = 0,        xend = mean(de), y = mean(dc),  yend = mean(dc)), linetype = 2, size = .25)
p

#' Panel 2: uncertainty cloud (PSA simulations shown as grey points).
p = ggplot(data = df, aes(de, dc)) +
  geom_point(col = "gray34", size = .4) +
  theme_classic() +
  theme(
    axis.line = element_blank(), axis.text.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks = element_blank(),
    legend.position = "none", panel.background = element_blank(),
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), plot.background = element_blank()
  ) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_segment(aes(x = 0,    y = -Inf, xend = 0,   yend = Inf), arrow = arrow(length = unit(0.30, "cm"), type = "closed")) +
  geom_segment(aes(x = -Inf, y = 0,    xend = Inf, yend = 0),   arrow = arrow(length = unit(0.30, "cm"), type = "closed")) +
  labs(title = "", y = "Cost differential", x = "Effectiveness differential") +
  theme(
    plot.title  = element_text(hjust = 0.5, size = 18),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15)
  ) +
  annotate("text", 0, Inf, label = "$\\Delta_c$", hjust = 2, vjust = 1, size = 5) +
  annotate("text", Inf, 0, label = "$\\Delta_e$", hjust = 1, vjust = 2, size = 5)
p

#' Panels 3 & 4: BCEA's ceplane.plot() with the sustainability line at two
#' different WTP thresholds (k=25000 and k=10000).
#' The sustainability area is the region below the line Delta_c = k * Delta_e;
#' PSA points in that area support the cost-effectiveness of the new intervention.
p = ceplane.plot(m, graph = "gg") +
  theme(
    axis.line = element_blank(), axis.text.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks = element_blank(),
    legend.position = "none", panel.background = element_blank(),
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), plot.background = element_blank()
  ) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_segment(aes(x = 0,    y = -Inf, xend = 0,   yend = Inf), arrow = arrow(length = unit(0.30, "cm"), type = "closed")) +
  geom_segment(aes(x = -Inf, y = 0,    xend = Inf, yend = 0),   arrow = arrow(length = unit(0.30, "cm"), type = "closed")) +
  geom_point(aes(x = mean(df$de), y = mean(df$dc)), col = "red", size = 2) +
  labs(title = "", y = "Cost differential", x = "Effectiveness differential") +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
        plot.title = element_text(face = "plain")) +
  annotate("text", 0, Inf, label = "$\\Delta_c$", hjust = 2, vjust = 1, size = 5) +
  annotate("text", Inf, 0, label = "$\\Delta_e$", hjust = 1, vjust = 2, size = 5)
p

p = ceplane.plot(m, graph = "gg", wtp = 10000) +
  theme(
    axis.line = element_blank(), axis.text.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks = element_blank(),
    legend.position = "none", panel.background = element_blank(),
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), plot.background = element_blank()
  ) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_segment(aes(x = 0,    y = -Inf, xend = 0,   yend = Inf), arrow = arrow(length = unit(0.30, "cm"), type = "closed")) +
  geom_segment(aes(x = -Inf, y = 0,    xend = Inf, yend = 0),   arrow = arrow(length = unit(0.30, "cm"), type = "closed")) +
  geom_point(aes(x = mean(df$de), y = mean(df$dc)), col = "red", size = 2) +
  labs(title = "", y = "Cost differential", x = "Effectiveness differential") +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
        plot.title = element_text(face = "plain")) +
  annotate("text", 0, Inf, label = "$\\Delta_c$", hjust = 2, vjust = 1, size = 5) +
  annotate("text", Inf, 0, label = "$\\Delta_e$", hjust = 1, vjust = 2, size = 5)
p


# ==============================================================================
# SECTION: Cost-Effectiveness Acceptability Curve (CEAC)
# ==============================================================================

#' The CEAC plots Pr(IB(theta) > 0) as a function of the willingness-to-pay k.
#' For each value of k, this is the proportion of PSA simulations that fall
#' inside the sustainability area of the CE plane -- i.e. where the new
#' intervention is cost-effective at that threshold.
#' BCEA computes and stores this automatically inside the bcea object.
ceac.plot(m, graph = "gg") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = "")


# ==============================================================================
# SECTION: Multi-intervention analysis (Smoking dataset)
# ==============================================================================

#' When T > 2, BCEA produces pairwise CEACs (one per comparator), the
#' probability that each intervention is the most cost-effective (via
#' multi.ce()), the Cost-Effectiveness Acceptability Frontier (CEAF) and the
#' Cost-Effectiveness Efficiency Frontier (CEEF).
#'
#' The Smoking dataset has four interventions; t=4 is set as the reference.
#' Kmax=1000 restricts the willingness-to-pay grid to [0, 1000].

data(Smoking)
m2 = bcea(
  eff, cost, ref = 4, Kmax = 1000,
  interventions = c("\\(t=1\\)", "\\(t=2\\)", "\\(t=3\\)", "\\(t=4\\)")
)

#' Pairwise CEACs: one curve per comparator, showing the probability that
#' t=4 is cost-effective against each alternative.
ceac.plot(m2, graph = "gg") +
  labs(title = "") + theme_bw() +
  theme(
    legend.position.inside = c(.75, .35), legend.position = "inside",
    legend.background = element_blank(), legend.text = element_text(size = 12)
  )

#' multi.ce() computes, for each value of k, the probability that each
#' intervention (not just the reference) is the most cost-effective overall.
mce = multi.ce(m2)

#' Plots the per-intervention probability of being most cost-effective.
#' The four curves sum to 1 at every value of k.
cc = ceac.plot(mce, graph = "gg")
cc$data |>
  mutate(comparison = paste0("\\(t=", comparison, "\\)") |> as.factor()) |>
  ggplot(aes(k, ceac, col = comparison)) +
  geom_line(linewidth = .95) +
  labs(title = "", color = "") +
  xlab("Willingness to pay") +
  ylab("Probability of most cost-effective") +
  scale_x_continuous(label = scales::comma, limits = c(0, 1035)) +
  directlabels::geom_dl(
    aes(label = comparison),
    method = list("last.points", cex = 1.2, colour = "black")
  ) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#000000", "#1F77B4", "#FF7F0E", "#2CA02C"))

#' CEAF: the probability of cost-effectiveness for the *optimal* intervention
#' at each value of k, tracing the upper envelope of the per-intervention curves.
#' "Switch points" where the optimal intervention changes appear as kinks.
ceaf.plot(mce, graph = "gg") +
  labs(title = "") + theme_bw() +
  geom_line(linewidth = .95) +
  scale_x_continuous(label = scales::comma)

#' CEEF: each intervention plotted as a point at its average cost and
#' effectiveness, with the efficiency frontier connecting the Pareto-optimal
#' interventions.  Interventions inside the frontier are never the most
#' cost-effective for any value of k.
#' dominance=FALSE suppresses dominated-intervention annotations.
ceef.plot(mce, graph = "gg", print.summary = FALSE, dominance = FALSE) +
  labs(title = "") +
  theme(legend.text = element_text(size = 12))
