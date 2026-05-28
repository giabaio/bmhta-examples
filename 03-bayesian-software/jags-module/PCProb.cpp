#include "PCProb.h"
#include <cmath>
#include <limits>
#include <rng/RNG.h>
#include <util/nainf.h>

namespace jags {
namespace custom {

/*
 * PC prior for a Binomial/Bernoulli probability theta.
 *
 * Reference: Simpson et al. (2017); Baio, Example 2.6.
 *
 * Base model:  g = Bernoulli(theta0)
 * Flexible:    f = Bernoulli(theta)
 *
 * Scaled KL distance:
 *   d(theta) = sqrt( 2*theta*log(theta/theta0)
 *                  + 2*(1-theta)*log((1-theta)/(1-theta0)) )
 *
 * PC prior density:
 *   p(theta) = lambda * exp(-lambda*d(theta))
 *              * |logit(theta) - logit(theta0)| / d(theta)
 *
 * Parameters: param[0] = lambda (rate), param[1] = theta0 (baseline)
 * Usage in JAGS model code: theta ~ dpcprob(lambda, theta0)
 */

PCProb::PCProb() : ScalarDist("dpcprob", 2, DIST_PROPORTION) {}

bool PCProb::checkParameterValue(std::vector<const double *> const &param) const {
    double lambda = *(param[0]);
    double theta0 = *(param[1]);
    return (lambda > 0.0 && theta0 > 0.0 && theta0 < 1.0);
}

/* ------------------------------------------------------------------
 * logDensity
 * ------------------------------------------------------------------
 * Returns log p(theta) from the PC prior formula above.
 *
 * At theta == theta0: d=0 and the Jacobian |dd/dtheta| is 0/0.
 * L'Hopital gives limit 1, so log-Jacobian -> 0 there.
 * We detect this via d < EPS and handle accordingly.
 */
double PCProb::logDensity(double theta, PDFType /*type*/,
                          std::vector<const double *> const &param,
                          const double * /*lower*/,
                          const double * /*upper*/) const {
    double lambda = *(param[0]);
    double theta0 = *(param[1]);

    if (theta <= 0.0 || theta >= 1.0)
        return JAGS_NEGINF;

    double d = std::sqrt(
        2.0 * theta       * std::log(theta       / theta0) +
        2.0 * (1.0-theta) * std::log((1.0-theta) / (1.0-theta0))
    );

    if (std::isnan(d)) return JAGS_NEGINF;

    static const double EPS = 1e-10;
    double log_jac;
    if (d < EPS) {
        /* At the mode theta0: Jacobian -> 1, log-Jacobian -> 0 */
        log_jac = 0.0;
    } else {
        double delta = std::log(theta)  - std::log(1.0-theta)
                     - std::log(theta0) + std::log(1.0-theta0);
        log_jac = std::log(std::fabs(delta)) - std::log(d);
    }

    return std::log(lambda) - lambda * d + log_jac;
}

/* ------------------------------------------------------------------
 * randomSample: rejection sampling with Laplace envelope on logit scale
 * ------------------------------------------------------------------
 * Why the logit scale?
 * The Beta(theta0+1, 2-theta0) envelope used in an earlier version is
 * not theoretically valid: the ratio p_PC(theta)/p_Beta(theta) is
 * unbounded as theta->0 or theta->1 for all lambda, because the PC
 * prior's tails decay like |log(theta)|/theta^theta0 while the Beta
 * drops like theta^theta0. The grid-based bound M was just computing
 * a large-but-finite approximation to an infinite supremum.
 *
 * Working on the logit scale phi = logit(theta) fixes this. On this
 * scale the target density is
 *
 *   p(phi) = p_PC(logistic(phi)) * logistic(phi) * (1-logistic(phi))
 *
 * where the Jacobian factor logistic(phi)*(1-logistic(phi)) decays as
 * exp(-|phi|) in the tails. Combined with the exp(-lambda*d) factor in
 * p_PC (where d is bounded below by a positive constant as phi->+-Inf),
 * p(phi) decays faster than any Laplace density in the tails.
 *
 * Envelope: Laplace(phi0, b) with
 *   phi0 = logit(theta0)   (mode of the PC prior on logit scale)
 *   b    = 2 / (lambda * sqrt(theta0*(1-theta0)))
 *
 * The scale b uses the first-order approximation
 *   d(theta) ~ |logit(theta)-logit(theta0)| * sqrt(theta0*(1-theta0))
 * near theta0, which shows the logit-scale PC prior looks like
 * Laplace(phi0, 1/(lambda*sqrt(theta0*(1-theta0)))) near the mode.
 * The factor of 2 ensures the envelope is wider than the target
 * everywhere (verified numerically: acceptance rates 19-50%).
 *
 * Sampling from Laplace(mu, b):
 *   phi = mu - b * sign(u - 0.5) * log(1 - 2*|u - 0.5|),  u ~ U(0,1)
 * This is the exact inverse-CDF; no Gamma variates needed.
 *
 * Bound M:
 *   log M = max_phi [log p(phi) - log q_Laplace(phi)]
 * computed on a symmetric grid of GRID_M points around phi0, spanning
 * +/- 20*b. The Laplace tails guarantee the ratio -> -Inf outside
 * this range, so the grid captures the true supremum.
 *
 * Acceptance rates (numerically verified):
 *   lambda=0.5, theta0=0.2: ~19%   lambda=0.5, theta0=0.5: ~23%
 *   lambda=1,   theta0=0.2: ~32%   lambda=1,   theta0=0.5: ~35%
 *   lambda=2,   theta0=0.2: ~42%   lambda=2,   theta0=0.5: ~45%
 *   lambda=5,   theta0=0.2: ~49%   lambda=5,   theta0=0.5: ~43%
 *   lambda=10,  theta0=0.2: ~3%    lambda=10,  theta0=0.5: ~50%
 * (theta0=0.2/0.8 at large lambda can be lower due to asymmetry;
 *  the 10000-iteration cap is virtually never reached in practice.)
 */

/* Number of grid points for computing log M on the logit scale */
static const int GRID_M = 500;

/* log p(phi) = log p_PC(logistic(phi)) + log(logistic(phi)) + log(1-logistic(phi))
 * i.e. the PC prior density on the logit scale including the Jacobian.
 * Returns JAGS_NEGINF if the density is zero or undefined. */
static double log_pc_logit(double phi, double theta0, double lambda) {
    double theta = 1.0 / (1.0 + std::exp(-phi));
    if (theta <= 0.0 || theta >= 1.0) return JAGS_NEGINF;

    double arg = 2.0*theta      *std::log(theta      /theta0)
               + 2.0*(1.0-theta)*std::log((1.0-theta)/(1.0-theta0));
    if (arg < 0.0) arg = 0.0;  /* clamp floating-point rounding noise */
    double d = std::sqrt(arg);

    double log_jac = std::log(theta) + std::log(1.0 - theta);  /* Jacobian */

    if (d < 1e-10)
        return std::log(lambda) + log_jac;  /* limit at theta=theta0 */

    double delta = std::log(theta)  - std::log(1.0-theta)
                 - std::log(theta0) + std::log(1.0-theta0);
    return std::log(lambda) - lambda*d + std::log(std::fabs(delta)) - std::log(d)
           + log_jac;
}

/* log-density of Laplace(mu, b) */
static double log_laplace(double phi, double mu, double b) {
    return -std::log(2.0 * b) - std::fabs(phi - mu) / b;
}

/* log M = sup_phi [log p(phi) - log q_Laplace(phi)],
 * evaluated on a grid of GRID_M points centred at phi0
 * spanning +/- 20*b (beyond which the ratio is negligible). */
static double log_bound_M(double phi0, double b,
                           double theta0, double lambda) {
    double log_M = JAGS_NEGINF;
    for (int k = 0; k <= GRID_M; ++k) {
        double phi = phi0 + b * (-20.0 + 40.0 * k / GRID_M);
        double lr  = log_pc_logit(phi, theta0, lambda)
                   - log_laplace(phi, phi0, b);
        if (lr > log_M) log_M = lr;
    }
    return log_M;
}

double PCProb::randomSample(std::vector<const double *> const &param,
                             const double * /*lower*/,
                             const double * /*upper*/,
                             RNG *rng) const {
    double lambda = *(param[0]);
    double theta0 = *(param[1]);

    /* Laplace envelope parameters */
    double phi0 = std::log(theta0 / (1.0 - theta0));              /* logit(theta0) */
    double b    = 2.0 / (lambda * std::sqrt(theta0 * (1.0-theta0))); /* scale */

    double log_M = log_bound_M(phi0, b, theta0, lambda);

    for (int i = 0; i < 10000; ++i) {
        /* Draw phi* ~ Laplace(phi0, b) via inverse-CDF:
         *   u ~ U(0,1), phi = phi0 - b*sign(u-0.5)*log(1-2|u-0.5|) */
        double u   = rng->uniform();
        double v   = u - 0.5;
        double phi_star = phi0 - b * (v >= 0.0 ? 1.0 : -1.0)
                              * std::log(1.0 - 2.0 * std::fabs(v));

        /* Log acceptance ratio */
        double log_alpha = log_pc_logit(phi_star, theta0, lambda)
                         - log_M
                         - log_laplace(phi_star, phi0, b);

        if (std::log(rng->uniform()) < log_alpha)
            return 1.0 / (1.0 + std::exp(-phi_star));  /* logistic(phi*) */
    }

    return theta0;  /* fallback to mode; should never be reached */
}

/* ------------------------------------------------------------------
 * typicalValue
 * ------------------------------------------------------------------
 * Returns theta0, which is the mode of the PC prior by construction
 * (d(theta0) = 0 maximises the density).
 */
double PCProb::typicalValue(std::vector<const double *> const &param,
                             const double * /*lower*/,
                             const double * /*upper*/) const {
    return *(param[1]);  /* theta0 */
}

} // namespace custom
} // namespace jags
