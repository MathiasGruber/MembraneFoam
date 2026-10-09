// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef membraneFluxModel_H
#define membraneFluxModel_H
#include <cmath>
#include <stdexcept>

namespace membrane
{
// SI units: A [m/(Pa s)], B [m/s], K [s/m], osmoticCoeff [Pa].
// AL-FS orientation, zero hydraulic pressure difference, draw >= feed.
// Multiplying the logarithmic model by K and exponentiating gives this
// residual without division by K or log(0). Its K=0 and B=0 limits are exact.
inline double residual(double j, double feed, double draw, double a,
                       double b, double k, double osmoticCoeff, bool advanced)
{
    const double salt = advanced ? b : 0.0;
    return j + salt + a*osmoticCoeff*feed
         - (salt + a*osmoticCoeff*draw)*std::exp(-j*k);
}

inline double flux(double feed, double draw, double a, double b, double k,
                   double osmoticCoeff, bool advanced, int* iterations=nullptr)
{
    if (!std::isfinite(feed+draw+a+b+k+osmoticCoeff) || feed < 0 || draw < feed
        || a <= 0 || b < 0 || k < 0 || osmoticCoeff <= 0)
        throw std::domain_error("FO flux requires finite A>0, B>=0, K>=0, osmoticCoeff>0 and 0<=feed<=draw");
    double lo = 0.0, hi = a*osmoticCoeff*(draw-feed);
    const double salt = advanced ? b : 0.0;
    const double feedTerm = salt + a*osmoticCoeff*feed;
    const double drawTerm = salt + a*osmoticCoeff*draw;
    double x = lo + 0.5*(hi-lo);
    double previousMagnitude = HUGE_VAL;
    int n = 0;
    // F is monotone and F' >= 1. Thus |F(x)| bounds the flux error.
    // Keep a bracket and bisect if Newton leaves it or stops reducing F.
    while (hi-lo > 1e-16 + 1e-12*hi && n < 100)
    {
        const double exponential = drawTerm*std::exp(-x*k);
        const double f = x + feedTerm - exponential;
        ++n;
        if (std::abs(f) <= 0.5*(1e-16 + 1e-12*x))
        {
            if (iterations) *iterations += n;
            return x;
        }
        if (f > 0) hi = x;
        else lo = x;
        const double candidate = x - f/(1 + k*exponential);
        x = candidate > lo && candidate < hi
            && std::abs(f) <= 0.5*previousMagnitude
            ? candidate : lo + 0.5*(hi-lo);
        previousMagnitude = std::abs(f);
    }
    if (iterations) *iterations += n;
    if (hi-lo > 1e-16 + 1e-12*hi)
        throw std::runtime_error("FO flux root did not converge within its bracket");
    return lo + 0.5*(hi-lo);
}
}
#endif
