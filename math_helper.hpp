#pragma once

#include <gsl/gsl_sf_gamma.h>

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

inline
REAL binomial_prefactor(int exponent, int pow1, int pow2, REAL x1, REAL x2)
{
    REAL s = 0;
    for(int i = 0; i < (1 + exponent); i++) {
        int j = exponent - i;
        if (i <= pow1 && j <= pow2) {
            s += binomial(pow1, i) * std::pow(x1, pow1 - i) * binomial(pow2, j) * std::pow(x2, pow2 - j);
        }
    }
    return s;
}

inline REAL 
boys(REAL n, REAL x) {
    if (x < 0.) {   throw;  }
    REAL radius = 0.0001;
    if (x < radius) {
        return 1. / (2.0 * n + 1.0);
    } else {
        REAL numerator = gsl_sf_gamma_inc_P(n+0.5, x) * gsl_sf_gamma(n+0.5);
        REAL denominator = 2.0 * std::pow(x, n+0.5);
        return numerator / denominator;
    }
}


