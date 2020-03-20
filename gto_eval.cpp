#include "common.hpp"
#include "gto.hpp"
#include "gto_eval.hpp"
#include "system.hpp"
#include <cmath>

#include <algorithm>
#include <boost/array.hpp>

namespace MOSolver {

static REAL 
_overlap1D(
        REAL exponent1, REAL x1, int l1, 
        REAL exponent2, REAL x2, int l2)
{
    REAL s(0);
    REAL gamma = exponent1 + exponent2;
    REAL center_x = (exponent1 * x1 + exponent2 * x2) / gamma;
    REAL Ra = center_x - x1;
    REAL Rb = center_x - x2;
    for(int i = 0; i < (1 + int(std::floor((l1+l2) / 2))); i++) {
        s += binomial_prefactor(2 * i, l1, l2, Ra, Rb) * factorial2(2*i-1) / std::pow(2.0*gamma, i);
    }
    return s;
}

static REAL 
_overlap3D(
        REAL exponent1, Vector3Real center1, int l1, int m1, int n1, 
        REAL exponent2, Vector3Real center2, int l2, int m2, int n2) {

    REAL dist2 = (center1 - center2).squaredNorm();
    REAL gamma = exponent1 + exponent2;
    if (gamma == REAL(0.)) {
        std::cout << __LINE__ << std::endl;
        throw;
    }
    REAL prefactor = std::pow(M_PI/gamma, 1.5) * std::exp(-exponent1 * exponent2 * dist2 / gamma);
    if (prefactor == REAL(0.)) {
        std::cout << __LINE__ << std::endl;
        throw;
    }

    REAL sx = _overlap1D(exponent1, center1[0], l1, exponent2, center2[0], l2);
    REAL sy = _overlap1D(exponent1, center1[1], m1, exponent2, center2[1], m2);
    REAL sz = _overlap1D(exponent1, center1[2], n1, exponent2, center2[2], n2);
    REAL ret = prefactor * sx * sy * sz;
    return ret;
}

REAL
overlap_PGTO(const PrimitiveGTO &lhs, const PrimitiveGTO &rhs) {
    REAL s = 0;
    s = _overlap3D(
            lhs.exponent, lhs.center, lhs.l, lhs.m, lhs.n,
            rhs.exponent, rhs.center, rhs.l, rhs.m, rhs.n );
    return s * lhs.norm_factor * rhs.norm_factor;
}

REAL
kinetic_PGTO(const PrimitiveGTO &lhs, const PrimitiveGTO &rhs) {
    REAL exp1 = lhs.exponent;
    REAL exp2 = rhs.exponent;
    Vector3Real Ra = lhs.center;
    Vector3Real Rb = rhs.center;

    int l1 = lhs.l; int m1 = lhs.m; int n1 = lhs.n;
    int l2 = rhs.l; int m2 = rhs.m; int n2 = rhs.n;

    REAL term1 = exp2 * (2*(l2+m2+n2)+3) * _overlap3D(exp1,Ra,l1,m1,n1,exp2,Rb,l2,m2,n2);

    REAL term2 = -2.0 * std::pow(exp2, 2) * (
            _overlap3D(exp1, Ra, l1, m1, n1, exp2, Rb, l2 + 2, m2, n2) +
            _overlap3D(exp1, Ra, l1, m1, n1, exp2, Rb, l2, m2 + 2, n2) + 
            _overlap3D(exp1, Ra, l1, m1, n1, exp2, Rb, l2, m2, n2 + 2) );

    REAL term3 = -0.5 * (
            l2*(l2-1) * _overlap3D(exp1, Ra, l1, m1, n1, exp2, Rb, l2-2, m2, n2) + 
            m2*(m2-1) * _overlap3D(exp1, Ra, l1, m1, n1, exp2, Rb, l2, m2-2, n2) + 
            n2*(n2-1) * _overlap3D(exp1, Ra, l1, m1, n1, exp2, Rb, l2, m2, n2-2) );

    return (term1 + term2 + term3) * lhs.norm_factor * rhs.norm_factor;
}

typedef boost::array<REAL, 2 * MAX_L + 1> g_list_array;

static void
_g_list(int l1, int l2, const REAL PA_x, const REAL PB_x, const REAL CP_x, REAL gamma, g_list_array &g_list_out)
{
    std::fill(g_list_out.begin(), g_list_out.end(), 0.);
    for(int i = 0; i < 1 + l1 + l2; i++) {
        for(int r = 0; r < int(1+std::floor(i/2.0)); r++) {
            for(int u = 0; u < int(1+std::floor((i-2*r)/2)); u++) {
                size_t I = i - 2*r - u;
                REAL term1 = std::pow(-1, i) * binomial_prefactor(i, l1, l2, PA_x, PB_x );
                REAL term2_numerator = std::pow(-1, u) * factorial(i) * std::pow(CP_x, i-2*r-2*u) * std::pow(0.25/gamma, r+u);
                REAL term2_denominator = factorial(r) * factorial(u) * factorial(i-2*r-2*u);
                g_list_out[I] = g_list_out[I] + term1 * term2_numerator / term2_denominator;
            }
        }
    }
}

REAL
nuclear_attraction_PGTO(const PrimitiveGTO &lhs, const PrimitiveGTO &rhs, const Atom &atom) {
    // <lhs | 1/Rzx | rhs>
    int l1 = lhs.l; int m1 = lhs.m; int n1 = lhs.n;
    int l2 = rhs.l; int m2 = rhs.m; int n2 = rhs.n;
    REAL exp1 = lhs.exponent;   
    REAL exp2 = rhs.exponent;
    REAL gamma = exp1 + exp2;
    REAL dist2 = (lhs.center - rhs.center).squaredNorm();
    Vector3Real Rp = (exp1 * lhs.center + exp2 * rhs.center) / gamma;
    REAL pc2 = (atom.center - Rp).squaredNorm();
    REAL pre = 2.0 * M_PI / gamma * std::exp(-exp1 * exp2 * dist2 / gamma);

    Vector3Real PA = Rp - lhs.center;
    Vector3Real PB = Rp - rhs.center;
    Vector3Real CP = Rp - atom.center;

    g_list_array a_x, a_y, a_z;
    _g_list(l1, l2, PA[0], PB[0], CP[0], gamma, a_x);
    _g_list(m1, m2, PA[1], PB[1], CP[1], gamma, a_y);
    _g_list(n1, n2, PA[2], PB[2], CP[2], gamma, a_z);
    REAL s = 0.;
    for(int I = 0; I < 1+l1+l2; I++) {
        for(int J = 0; J < 1+m1+m2; J++) {
            for(int K = 0; K < 1+n1+n2; K++) {
                s += a_x[I] * a_y[J] * a_z[K] * boys(I+J+K, gamma*pc2);
            }
        }
    }
    s = s * pre * lhs.norm_factor * rhs.norm_factor;
    return s;
}

typedef boost::array<REAL, 4*MAX_L + 1> c_list_array;

// Returns the C[I] in the paper pp.2320.
static void
_c_list(REAL exp1, REAL x1, int l1, REAL exp2, REAL x2, int l2, REAL px,
    REAL exp3, REAL x3, int l3, REAL exp4, REAL x4, int l4, REAL qx, c_list_array &c_list_out)
{
    std::fill(c_list_out.begin(), c_list_out.end(), 0.);
    //std::cout << " l1: " << l1 << " l2: " << l2 << " l3: " << l3 << " l4: " << l4 << std::endl;
    REAL pa = px - x1;
    REAL pb = px - x2;
    REAL qc = qx - x3;
    REAL qd = qx - x4;
    REAL gamma1 = exp1 + exp2;
    REAL gamma2 = exp3 + exp4;
    REAL delta = 0.25*(1./gamma1 + 1./gamma2);
    for(int i1 = 0; i1 < l1 + l2 + 1; i1++) {
        for(int i2 = 0; i2 < l3 + l4 + 1; i2++) {
            for(int r1 = 0; r1 < int(std::floor(i1/2.)) + 1; r1++) {
                for(int r2 = 0; r2 < int(std::floor(i2/2.)) + 1; r2++) {
                    for(int u = 0; u < std::floor((i1+i2)/2-r1-r2 + 1); u++) {
                        REAL f1 = binomial_prefactor(i1,l1,l2,pa,pb) * factorial(i1) * std::pow(4*gamma1,r1) / (std::pow(4*gamma1,i1) * factorial(r1) * factorial(i1-2*r1));
                        REAL f2 = binomial_prefactor(i2,l3,l4,qc,qd) * factorial(i2) * std::pow(4*gamma2,r2) / (std::pow(4*gamma2,i2) * factorial(r2) * factorial(i2-2*r2));
                        REAL f3_numerator = factorial(i1+i2-2*(r1+r2)) * std::pow(-1., u) * std::pow(qx-px, i1+i2-2*(r1+r2)-2*u);
                        REAL f3_denominator=factorial(u) * factorial(i1+i2-2*(r1+r2)-2*u) * std::pow(delta, i1+i2-2*(r1+r2)-u);
                        int I = i1+i2-2*(r1+r2)-u;
                        //c_list[I] += f1*f2*f3_numerator/f3_denominator;
                        c_list_out[I] += f1*std::pow(-1., i2)*f2*f3_numerator/f3_denominator;
                    }
                }
            }
        }
    }
    return;
}

REAL
electron_reulsion_PGTO(const PrimitiveGTO& pgto1, const PrimitiveGTO &pgto2, const PrimitiveGTO &pgto3, const PrimitiveGTO &pgto4)
{
    REAL gamma1 = pgto1.exponent + pgto2.exponent;
    REAL gamma2 = pgto3.exponent + pgto4.exponent;

    Vector3Real p = (pgto1.exponent * pgto1.center + pgto2.exponent * pgto2.center) / gamma1;
    Vector3Real q = (pgto3.exponent * pgto3.center + pgto4.exponent * pgto4.center) / gamma2;
    REAL rab2 = (pgto1.center - pgto2.center).squaredNorm();
    REAL rcd2 = (pgto3.center - pgto4.center).squaredNorm();
    REAL rpq2 = (p - q).squaredNorm();

    REAL delta = (1./gamma1 + 1./gamma2) / 4.;
    REAL s(0.);
    c_list_array c_list_x, c_list_y, c_list_z;
    _c_list(pgto1.exponent, pgto1.center[0], pgto1.l, pgto2.exponent, pgto2.center[0], pgto2.l, p[0], 
            pgto3.exponent, pgto3.center[0], pgto3.l, pgto4.exponent, pgto4.center[0], pgto4.l, q[0], c_list_x );
    _c_list(pgto1.exponent, pgto1.center[1], pgto1.m, pgto2.exponent, pgto2.center[1], pgto2.m, p[1], 
            pgto3.exponent, pgto3.center[1], pgto3.m, pgto4.exponent, pgto4.center[1], pgto4.m, q[1], c_list_y );
    _c_list(pgto1.exponent, pgto1.center[2], pgto1.n, pgto2.exponent, pgto2.center[2], pgto2.n, p[2], 
            pgto3.exponent, pgto3.center[2], pgto3.n, pgto4.exponent, pgto4.center[2], pgto4.n, q[2], c_list_z );

    for(int I = 0; I < 1 + pgto1.l + pgto2.l + pgto3.l + pgto4.l; I++) {
        for(int J = 0; J < 1 + pgto1.m + pgto2.m + pgto3.m + pgto4.m; J++) {
            for(int K = 0; K < 1 + pgto1.n + pgto2.n + pgto3.n + pgto4.n; K++) {
                s += c_list_x[I] * c_list_y[J] * c_list_z[K] * boys(I+J+K, rpq2/(4.*delta));
            }
        }
    }
    s *= 2.0*std::pow(M_PI, 2.5)/ (gamma1 * gamma2)/ std::pow(gamma1 + gamma2, 0.5) * \
        std::exp(-(pgto1.exponent*pgto2.exponent*rab2/gamma1) - (pgto3.exponent*pgto4.exponent*rcd2/gamma2));

    s *= pgto1.norm_factor * pgto2.norm_factor * pgto3.norm_factor * pgto4.norm_factor;
    return s;
}

REAL
overlap_CGTO(const ContractedGTO &lhs, const ContractedGTO &rhs) {
    REAL s(0.);
    size_t npgto_lhs(lhs.num_pgtos());
    size_t npgto_rhs(rhs.num_pgtos());
    for(size_t i = 0; i < npgto_lhs; i++) {
        for(size_t j = 0; j < npgto_rhs; j++) {
            s += lhs.coeff_list[i] * rhs.coeff_list[j] * overlap_PGTO(lhs.pgto_list[i], rhs.pgto_list[j]);
        }
    }
    s = s * lhs.norm_factor * rhs.norm_factor;
    return s;
}

REAL
kinetic_CGTO(const ContractedGTO &lhs, const ContractedGTO &rhs) {
    REAL t(0.);
    size_t npgto_lhs(lhs.num_pgtos());
    size_t npgto_rhs(rhs.num_pgtos());
    for(size_t i = 0; i < npgto_lhs; i++) {
        for(size_t j = 0; j < npgto_rhs; j++) {
            t += lhs.coeff_list[i] * rhs.coeff_list[j] * kinetic_PGTO(lhs.pgto_list[i], rhs.pgto_list[j]);
        }
    }
    t = t * lhs.norm_factor * rhs.norm_factor;
    return t;
}

REAL
nuclear_attraction_CGTO(const ContractedGTO &lhs, const ContractedGTO &rhs, const Atom &atom) {
    // <lhs | Z/Rzx | rhs>
    REAL t(0.);
    for (size_t i = 0; i < lhs.num_pgtos(); i++) {
        for(size_t j = 0; j < rhs.num_pgtos(); j++) {
            t += lhs.coeff_list[i] * rhs.coeff_list[j] * nuclear_attraction_PGTO(lhs.pgto_list[i], rhs.pgto_list[j], atom);
        }
    }
    t = -atom.atomic_number * t * lhs.norm_factor * rhs.norm_factor;
    return t;
}

REAL 
electron_repulsion_CGTO(const ContractedGTO &cgto1, const ContractedGTO &cgto2, const ContractedGTO &cgto3, const ContractedGTO &cgto4) 
{
    // < cgto1 cgto2 | 1/r | cgto3 cgto4 >
    REAL t(0.);
    for (size_t i = 0; i < cgto1.num_pgtos(); i++) {
        for(size_t j = 0; j < cgto2.num_pgtos(); j++) {
            for(size_t k = 0; k < cgto3.num_pgtos(); k++) {
                for(size_t l = 0; l < cgto4.num_pgtos(); l++) {
                    REAL e = electron_reulsion_PGTO(cgto1.pgto_list[i], cgto2.pgto_list[j], cgto3.pgto_list[k], cgto4.pgto_list[l] );
                    REAL prod_coeff = cgto1.coeff_list[i] * cgto2.coeff_list[j] * cgto3.coeff_list[k] * cgto4.coeff_list[l];
                    t += prod_coeff * e;
                }
            }
        }
    }
    return t * cgto1.norm_factor * cgto2.norm_factor * cgto3.norm_factor * cgto4.norm_factor;
}

} // namespace MOSolver 
