/*
 * PCProb.h
 * --------
 * Header for the PCProb class: PC prior for a Binomial/Bernoulli probability
 * parameter theta, implemented as a custom JAGS distribution.
 *
 * Base class: ScalarDist  (distribution/ScalarDist.h)
 * This is the correct base for distributions that implement their own
 * logDensity/randomSample/typicalValue rather than delegating to libRmath.
 * (RScalarDist, by contrast, is for distributions built on top of R's
 * d/p/q/r functions and has a different interface.)
 *
 * Pure virtuals that must be implemented:
 *   logDensity()          -- log p(x | params)
 *   randomSample()        -- draw one sample (used for initial values)
 *   typicalValue()        -- a high-density point (mode, mean, ...)
 *   checkParameterValue() -- validate parameter ranges
 *
 * ScalarDist already provides isSupportFixed() for DIST_PROPORTION,
 * so we do not need to implement it.
 *
 * Parameters (in order as passed by JAGS):
 *   param[0] = lambda   -- rate / penalisation factor  (lambda > 0)
 *   param[1] = theta0   -- baseline probability        (0 < theta0 < 1)
 *
 * Usage in JAGS model code after load.module("pcprob"):
 *   theta ~ dpcprob(lambda, theta0)
 */

#ifndef PCPROB_H_
#define PCPROB_H_

#include <distribution/ScalarDist.h>

namespace jags {
namespace custom {

class PCProb : public ScalarDist {
public:
    /* Constructor: registers the distribution as "dpcprob" with 2 parameters
     * and support DIST_PROPORTION (open interval (0,1)). */
    PCProb();

    /* logDensity()
     * ------------
     * Returns log p(theta | lambda, theta0) under the PC prior, or the
     * density itself depending on the PDFType argument (see below).
     * Returns -Inf (or 0 in density scale) outside (0,1).
     *
     * PDFType is an enum defined in distribution/Distribution.h:
     *   PDF_FULL        -- full log-density including all constants
     *   PDF_PRIOR       -- may omit terms depending only on parameters
     *   PDF_LIKELIHOOD  -- may omit terms depending only on x
     * We always compute the full density and ignore the type argument. */
    double logDensity(double x, PDFType type,
                      std::vector<double const *> const &parameters,
                      double const *lower, double const *upper) const override;

    /* randomSample()
     * ---------------
     * Draws one exact sample from the PC prior using rejection sampling
     * with a Beta(theta0+1, 2-theta0) envelope.  Called by JAGS to
     * generate initial values; not on the MCMC hot path. */
    double randomSample(std::vector<double const *> const &parameters,
                        double const *lower, double const *upper,
                        RNG *rng) const override;

    /* typicalValue()
     * ---------------
     * Returns theta0, the mode of the PC prior (d(theta0)=0 => maximum
     * density by construction). Used by JAGS as a starting-point hint. */
    double typicalValue(std::vector<double const *> const &parameters,
                        double const *lower, double const *upper) const override;

    /* checkParameterValue()
     * ----------------------
     * Returns true iff lambda > 0 and 0 < theta0 < 1. */
    bool checkParameterValue(
        std::vector<double const *> const &parameters) const override;
};

} // namespace custom
} // namespace jags

#endif /* PCPROB_H_ */
