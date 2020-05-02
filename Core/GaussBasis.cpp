#include "Core.hpp"
#include "GaussBasis.hpp"
#include <cmath>
#include "math_helper.hpp"

namespace MOSolver2 {

static REAL 
self_overlap1D(const REAL exp1, const REAL exp2, const int l)
{
    REAL s(0.);
    REAL gamma = exp1 + exp2;
    for(int i = 0; i < (1 + int(std::floor(l))); i++) {
        s += binomial_prefactor(2*i,l,l,0.,0.) * factorial2(2*i-1)/std::pow(2.0*gamma,i);
    }
    return s;
}

static REAL 
self_overlap(const REAL exp1, const REAL exp2, const int l, const int m, const int n)
{
    REAL gamma = exp1 + exp2;
    REAL prefactor = std::pow(M_PI/gamma, 1.5);
    REAL sx = self_overlap1D(exp1, exp2, l);
    REAL sy = self_overlap1D(exp1, exp2, m);
    REAL sz = self_overlap1D(exp1, exp2, n);
    REAL ret = prefactor * sx * sy * sz;
    return ret;
}

void
GaussBasis::normalize()
{
    // This function is basically the same as overlap defined in libGaussEval.
    for(size_t i = 0; i < this->length(); i++) {
        REAL v2 = self_overlap(exponents_[i], exponents_[i], this->l, this->m, this->n);
        this->primitive_normalization_factor_[i] = 1./std::sqrt(v2);
    }
    REAL total_norm2 = 0.;
    for(size_t i = 0; i < this->length(); i++) {
        for(size_t j = 0; j < this->length(); j++) {
            REAL c_i = coefficients_[i];
            REAL c_j = coefficients_[j];
            REAL prim_norm_i = primitive_normalization_factor_[i];
            REAL prim_norm_j = primitive_normalization_factor_[j];
            REAL exp_i = exponents_[i];
            REAL exp_j = exponents_[j];
            total_norm2 += prim_norm_i * prim_norm_j * c_i * c_j * self_overlap(exp_i, exp_j, this->l, this->m, this->n);
        }
    }
    this->normalization_factor_ = 1./std::sqrt(total_norm2);
}

}


