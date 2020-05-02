#pragma once

#include <vector>
#include "Definition.hpp"

namespace MOSolver2 {

struct GaussBasis {
    GaussBasis(
            const std::vector<REAL> &exponents, 
            const std::vector<REAL> &coefficients,
            int l, int m, int n, 
            const REAL x, const REAL y, const REAL z):
        exponents_(exponents), coefficients_(coefficients), 
        l(l), m(m), n(n), x(x), y(y), z(z)
    {
        if (exponents_.size() != coefficients_.size()) { throw; }
        this->primitive_normalization_factor_.resize(this->length());
        normalize();
    }

    void normalize();
    REAL value(const REAL x, const REAL y, const REAL z);
    size_t length() const {return this->exponents_.size();}

    //members variables;
    std::vector<REAL> exponents_;
    std::vector<REAL> coefficients_;
    std::vector<REAL> primitive_normalization_factor_;
    int  l, m, n;
    REAL x, y, z;
    REAL normalization_factor_;
};

struct BasisFunctions {
    void push_back(const GaussBasis &gf) {
        this->basis_function_list_.push_back(gf);
    }
    size_t size(void) const { return this->basis_function_list_.size();}
    const GaussBasis &operator[](const size_t index) const {
        return this->basis_function_list_[index];
    }
    GaussBasis &operator[](const size_t index) {
        return this->basis_function_list_[index];
    }
    std::vector<GaussBasis> basis_function_list_;
};

} // namespace
