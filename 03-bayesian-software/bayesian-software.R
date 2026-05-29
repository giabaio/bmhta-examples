#' ---
#' title: "Bayesian Software"
#' desc:  "Covers running JAGS from R for Bayesian inference, post-processing
#'         MCMC output (traceplots, autocorrelation, convergence diagnostics),
#'         key differences between BUGS and JAGS defaults, and using the
#'         zero-trick to implement non-standard distributions in JAGS."
#' ---


# ==============================================================================
# SECTION: Package installation
# ==============================================================================

#' R interfaces to JAGS and OpenBUGS respectively.
#' Only one is strictly needed; JAGS is the primary tool used in this book.
install.packages("R2jags")
install.packages("R2OpenBUGS")


# ==============================================================================
# SECTION: Drug example -- running a conjugate model in JAGS
# ==============================================================================

#' The model is a Beta-Binomial with a predictive node.
#' We observe y=15 successes out of m=20 trials and want to:
#'   (a) update the prior Beta(9.2, 13.8) to the posterior for theta, and
#'   (b) propagate that posterior uncertainty into a predictive distribution
#'       for y.pred, the number of successes in a future trial of size n=40.
#' P.crit flags whether y.pred meets or exceeds a critical threshold (ncrit=25).

# JAGS model code (save as "drug-MCMC.txt" or pass as an R function -- see below)
#
# model {
#   theta   ~ dbeta(a, b)                  # prior for the success probability
#   y       ~ dbin(theta, m)               # likelihood: observed successes
#   y.pred  ~ dbin(theta, n)               # posterior predictive distribution
#   P.crit  <- step(y.pred - ncrit + 0.5)  # 1 if y.pred >= ncrit, 0 otherwise
# }
#
# Note: BUGS/JAGS parameterises Normal distributions by precision (= 1/variance),
# not standard deviation.  Assignments use '<-', not '='.

#' Passing data to JAGS as a named list keeps model code general and values
#' easy to change without touching the model specification.
data = list(
  a     = 9.2,  # Beta prior: first shape parameter
  b     = 13.8, # Beta prior: second shape parameter
  y     = 15,   # observed successes
  m     = 20,   # observed trial size
  n     = 40,   # future trial size for the predictive distribution
  ncrit = 25    # critical threshold for P.crit
)

#' Alternatively, model code can be embedded directly as an R function
#' (avoids managing a separate .txt file).
#' R2jags converts it to a temporary file via R2WinBUGS::write.model().
model.code = function() {
  theta   ~ dbeta(a, b)
  y       ~ dbin(theta, m)
  y.pred  ~ dbin(theta, n)
  P.crit  <- step(y.pred - ncrit + 0.5)
}

#' Runs JAGS in the background and returns an R object with the MCMC output.
#' Key arguments:
#'   data              -- named list of data and constants
#'   parameters.to.save-- nodes to monitor (stored for post-processing)
#'   inits             -- NULL lets JAGS choose starting values automatically
#'   model.file        -- path to .txt file, or an R function (as here)
#'   n.chains          -- number of parallel Markov chains
#'   n.iter            -- total iterations per chain (includes burn-in)
#'   n.burnin          -- iterations discarded before sampling begins
#'   n.thin            -- retain one in every n.thin iterations (reduces storage)
#'   DIC               -- compute Deviance Information Criterion if TRUE
library(R2jags)
model = jags(
  data               = data,
  parameters.to.save = c("y.pred", "theta", "P.crit"),
  inits              = NULL,
  model.file         = model.code,
  n.chains           = 2,
  n.iter             = 6000,
  n.burnin           = 1000,
  n.thin             = 1,
  DIC                = TRUE
)


# ==============================================================================
# SECTION: Post-processing JAGS output
# ==============================================================================

#' Prints a summary table of posterior means, standard deviations, quantiles,
#' R-hat and effective sample size for all monitored nodes.
print(model, digit = 3)

#' Traceplots for selected nodes.
#' Good mixing shows the two chains overlapping immediately and staying
#' interleaved throughout -- here expected because the model is conjugate.
bmhe::traceplot(model, "theta")
bmhe::traceplot(model, "y.pred")

#' Autocorrelation function for monitored nodes.
#' Low autocorrelation at small lags indicates the chain mixes well and
#' successive draws are approximately independent.
bmhe::acfplot(model, parameter = c("theta", "y.pred"), col = "red")

#' Visual summary of the two main formal convergence diagnostics:
#'   R-hat (potential scale reduction): values close to 1 indicate convergence.
#'   n.eff (effective sample size):     values close to n.sims indicate low autocorrelation.
bmhe::diagplot(model, what = "Rhat")
bmhe::diagplot(model, what = "n.eff")

#' Posterior density plots for theta (continuous, in [0,1]) and
#' y.pred (discrete count in [0, 40]).  plot="bar" produces histograms.
bmhe::posteriorplot(model, plot = "bar", parameter = "theta")
bmhe::posteriorplot(model, plot = "bar", parameter = "y.pred")


# ==============================================================================
# SECTION: Differences between BUGS and JAGS
# ==============================================================================

#' Illustrates how the default settings of R2OpenBUGS::bugs() and R2jags::jags()
#' differ, even when called with nominally identical arguments.
#' Uses the built-in "eight schools" hierarchical model from R2OpenBUGS.

library(R2OpenBUGS)
library(R2jags)
set.seed(1)

#' Default function signatures -- shown for reference only (eval=FALSE).
#' Key differences to note:
#'   n.burnin: both default to floor(n.iter/2)
#'   n.thin:   BUGS defaults to 1; JAGS defaults to max(1, floor((n.iter-n.burnin)/1000))
#'   seed:     BUGS uses 1; JAGS uses 123
R2OpenBUGS::bugs(
  data, inits, parameters.to.save, n.iter, model.file = "model.txt",
  n.chains = 3, n.burnin = floor(n.iter / 2), n.thin = 1, saveExec = FALSE,
  restart = FALSE, debug = FALSE, DIC = TRUE, digits = 5, codaPkg = FALSE,
  OpenBUGS.pgm = NULL, working.directory = NULL, clearWD = FALSE,
  useWINE = FALSE, WINE = NULL, newWINE = TRUE, WINEPATH = NULL,
  bugs.seed = 1, summary.only = FALSE,
  save.history = (.Platform$OS.type == "windows" | useWINE == TRUE),
  over.relax = FALSE
)

R2jags::jags(
  data, inits, parameters.to.save, model.file = "model.bug",
  n.chains = 3, n.iter = 2000, n.burnin = floor(n.iter / 2),
  n.thin = max(1, floor((n.iter - n.burnin) / 1000)), DIC = TRUE,
  working.directory = NULL, jags.seed = 123, refresh = n.iter / 50,
  progress.bar = "text", digits = 5,
  RNGname = c(
    "Wichmann-Hill", "Marsaglia-Multicarry",
    "Super-Duper", "Mersenne-Twister"
  ),
  jags.module = c("glm", "dic"), quiet = FALSE
)

#' Loads the eight-schools model bundled with R2OpenBUGS and sets up the data.
model.file = system.file(package = "R2OpenBUGS", "model", "schools.txt")
data(schools)
J       = nrow(schools)
y       = schools$estimate
sigma.y = schools$sd
data    = list("J", "y", "sigma.y")

#' Initial value function: draws random starting points for each chain.
inits = function() {
  list(
    theta       = rnorm(J, 0, 100),
    mu.theta    = rnorm(1, 0, 100),
    sigma.theta = runif(1, 0, 100)
  )
}
parameters = c("theta", "mu.theta", "sigma.theta")

#' Runs both BUGS and JAGS with identical top-level arguments.
#' Despite n.iter=5000 and n.chains=3 in both cases, the retained n.sims differs:
#'   BUGS: n.thin=1 by default -> retains all (n.iter - n.burnin) * n.chains = 7500 sims
#'   JAGS: n.thin=max(1, floor((5000-2500)/1000))=2 by default -> retains 3750 sims
mbugs = bugs(data, inits, parameters, model.file, n.chains = 3, n.iter = 5000)
mjags = jags(data, inits, parameters, model.file, n.chains = 3, n.iter = 5000)

#' R2jags wraps all MCMC results inside a nested $BUGSoutput slot;
#' R2OpenBUGS stores results directly at the top level of the object.
#' bmhe helper functions handle both transparently.
names(mbugs)
names(mjags)

#' Accessing simulations from the JAGS object requires the extra level:
mjags$BUGSoutput$sims.matrix |> head() |> round(4)
# Equivalent for BUGS: mbugs$sims.matrix |> head()

#' Summary tables show the n.sims difference caused by the default thinning.
print(mbugs, digits = 3)
print(mjags, digits = 3)


# ==============================================================================
# SECTION: Zero-trick -- implementing a PC prior in JAGS
# ==============================================================================

#' JAGS/BUGS do not natively support the PC prior distribution derived in the
#' bayesian-computation chapter.  The "zero-trick" works around this:
#'
#'   1. Introduce a fake observation w = 0, modelled as w ~ Poisson(phi).
#'   2. Set phi = -log(target_density) + C, where C is a large positive constant
#'      that ensures phi > 0 (required because phi is a Poisson mean).
#'   3. The Poisson likelihood contribution exp(-phi) then equals the target
#'      density (up to the constant C), effectively imposing the desired prior.
#'
#' Here the target density is the PC prior for a Bernoulli parameter theta,
#' with baseline probability theta0 = 0.5 and scaling factor lambda = 1.

model_pc_bern = function() {
  # Pseudo-observation w=0; its Poisson likelihood encodes the PC prior
  w ~ dpois(phi)

  # KLD-based distance from the base model (theta0) and its derivative
  d     <- pow(
    2 * theta * log(theta / theta0) + 2 * (1 - theta) * log((1 - theta) / (1 - theta0)),
    0.5
  )
  deriv <- abs((log(theta) - log(1 - theta) - log(theta0) + log(1 - theta0)) / d)

  # phi must be positive: C=10000 is a large constant that acts as a scaling shift
  phi <- -log(lambda * exp(-lambda * d) * deriv) + C

  # Flat Beta(1,1) prior for theta so that phi carries all the likelihood weight
  theta ~ dbeta(1, 1)

  # Optional: y can be included if observed data are available;
  # set y=NA to run the model with no data (prior check only)
  y ~ dbin(theta, 10)
}

#' Runs the model with no observed data (y=NA) to verify the implied prior.
#' Two chains are started from opposite ends of [0,1] to check mixing.
#' n.iter=100000 and n.thin=4 are used because the zero-trick can produce
#' slow chain mixing, so a long run is advisable.
test = R2jags::jags(
  data               = list(y = NA, w = 0, C = 10000, lambda = 1, theta0 = 0.5),
  parameters.to.save = c("theta"),
  model.file         = model_pc_bern,
  n.chains           = 2,
  n.thin             = 4,
  n.iter             = 100000,
  DIC                = FALSE,
  inits              = list(list(theta = 0.1), list(theta = 0.9))
)

#' Summary of the posterior (here = prior, since y is unobserved).
print(test, digits = 3)

#' Monte Carlo estimate of Pr(theta > 0.75) under the PC prior.
#' Should agree with the numerical integration result (~0.22) from the
#' bayesian-computation chapter.
sum(test$BUGSoutput$sims.list$theta > 0.75) / test$BUGSoutput$n.sims

#' Plots the empirical density of the simulated PC prior as a bar chart.
tibble(x = test$BUGSoutput$sims.list$theta) |>
  ggplot(aes(x)) + geom_bar() +
  scale_x_binned(limits = c(0, 1)) +
  xlab("$\\theta$") + ylab("PC prior")
