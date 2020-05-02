#pragma once

#include "libGaussEval.hpp"
#include <boost/math/special_functions/gamma.hpp>
// #include <gsl/gsl_sf_gamma.h>  

namespace libGaussEval {

inline
REAL squaredNorm(
        const REAL x1, const REAL y1, const REAL z1,
        const REAL x2, const REAL y2, const REAL z2 )
{
    const REAL dx = x2 - x1;
    const REAL dy = y2 - y1;
    const REAL dz = z2 - z1;
    return std::pow(dx,2) + std::pow(dy,2) + std::pow(dz,2);
}

inline
int factorial(int n)
{
    if (n <= 1) {
        return 1;
    } else {
        return n * factorial(n-1);
    }
}

inline
int factorial2(int n)
{
    if (n <= 1) {   
        return 1;   
    } else {   
        return n * factorial2(n-2); 
    }
}

inline
int binomial(int n, int m) {
    // return the nCm
    return factorial(n) / factorial(m) / factorial(n-m);
}

inline REAL
binomial_prefactor(int exponent, int pow1, int pow2, REAL x1, REAL x2)
{
    REAL sum = 0.;
    for(int i = 0; i <= exponent; i++) {
        int j = exponent - i;
        if (i <= pow1 && j <= pow2) {
            sum += binomial(pow1, i) * binomial(pow2, j) * std::pow(x1, pow1-i) * std::pow(x2, pow2-j);
        }
    }
    return sum;
}

inline REAL 
boys(REAL n, REAL x) {
    if (x < 0.) {   throw;  }
    REAL radius = 0.0001;
    if (x < radius) {
        return 1. / (2.0 * n + 1.0);
    } else {
        // XXX At first the GSL was used for the incomplete gamma.
        // REAL numerator = gsl_sf_gamma_inc_P(n+0.5, x) * gsl_sf_gamma(n+0.5);
        REAL numerator = boost::math::gamma_p(n+0.5, x) * boost::math::tgamma(n+0.5);
        REAL denominator = 2.0 * std::pow(x, n+0.5);
        return numerator / denominator;
    }
}

}

