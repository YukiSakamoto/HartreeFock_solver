#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include "common.hpp"
#include "math_helper.hpp"

#include <iostream>

// defined in the gto_eval.hpp
struct ContractedGTO;
REAL overlap_CGTO(const ContractedGTO &lhs, const ContractedGTO &rhs); // gto_eval

struct PrimitiveGTO {
    PrimitiveGTO(REAL exponent, 
            int l, int m, int n, 
            REAL x, REAL y, REAL z):
        exponent(exponent), l(l), m(m), n(n), x(x), y(y), z(z), norm_factor(REAL(0)), center(x,y,z)
    {
        normalize();
    }

    REAL value(REAL rx, REAL ry, REAL rz);
    REAL value(Vector3Real v);
    void normalize();

    int l, m, n;
    REAL exponent;
    REAL x, y, z;   // the cental coordinate of the gaussian.
    Vector3Real center;
    REAL  norm_factor;
};


struct ContractedGTO {
    typedef std::vector<PrimitiveGTO>  PGTO_list;
    typedef std::vector<REAL>  Coeff_list;

    ContractedGTO(
            int l, int m, int n, 
            REAL x, REAL y, REAL z):
        l(l), m(m), n(n), x(x), y(y), z(z), norm_factor(REAL(0)), center(x,y,z), normalized(false)
    {}

    void add_primitiveGTO(REAL coeff, REAL exponent) {
        PrimitiveGTO pgto_new(exponent, l, m, n, x, y, z);
        pgto_list.push_back(pgto_new);
        coeff_list.push_back(coeff);
        normalized = false;
    }

    size_t num_pgtos(void) const
    {
        return pgto_list.size();
    }

    void normalize(void) {
        normalized = false;
        norm_factor = 1.;

        REAL s(overlap_CGTO(*this, *this));
        if (s == REAL(0.)) {
            std::cout << __LINE__ << std::endl;
            throw;
        }
        REAL norm( 1. / std::sqrt(s) );
        norm_factor = norm;
        normalized = true;
    }

    REAL value(REAL x, REAL y, REAL z) {
        Vector3Real r = Vector3Real(x,y,z);
        return value(r);
    }
    REAL value(Vector3Real r) {
        REAL acc = 0.;
        for(int i = 0; i < num_pgtos(); i++) {
            acc += pgto_list[i].value(r) * coeff_list[i];
        }
        return acc * norm_factor;
    }
    
    int l, m, n;
    REAL x, y, z;   // the cental coordinate of the gaussin.
    Vector3Real center;
    REAL  norm_factor;
    PGTO_list pgto_list;
    Coeff_list coeff_list;
    bool normalized;
};
