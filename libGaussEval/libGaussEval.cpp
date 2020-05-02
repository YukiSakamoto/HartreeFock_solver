#include <array>
#include <cmath>
#include "libGaussEval.hpp"
#include "libGaussUtilFunc.hpp"

namespace libGaussEval {

static REAL 
_overlap1D(
        const REAL exp1, const REAL x1, const int l1, 
        const REAL exp2, const REAL x2, const int l2)
{
    REAL s(0);
    REAL gamma = exp1 + exp2;
    REAL center_x = (exp1 * x1 + exp2 * x2) / gamma;
    REAL Ra = center_x - x1;
    REAL Rb = center_x - x2;
    for(int i = 0; i < (1 + int(std::floor((l1+l2) / 2))); i++) {
        s += binomial_prefactor(2 * i, l1, l2, Ra, Rb) * factorial2(2*i-1) / std::pow(2.0*gamma, i);
    }
    return s;
}

static REAL
_overlap3D(
        const REAL exp1, const int l1, const int m1, const int n1, 
        const REAL x1, const REAL y1, const REAL z1, 
        const REAL exp2, const int l2, const int m2, const int n2, 
        const REAL x2, const REAL y2, const REAL z2) 
{
    REAL dist2 = squaredNorm(x1, y1, z1, x2, y2, z2);
    REAL gamma = exp1+ exp2;
    if (gamma == REAL(0.)) {
        std::cout << __LINE__ << std::endl;
        throw;
    }
    REAL prefactor = std::pow(M_PI/gamma, 1.5) * std::exp(-exp1 * exp2 * dist2 / gamma);
    if (prefactor == REAL(0.)) {
        std::cout << __LINE__ << std::endl;
        throw;
    }

    REAL sx = _overlap1D(exp1, x1, l1, exp2, x2, l2);
    REAL sy = _overlap1D(exp1, y1, m1, exp2, y2, m2);
    REAL sz = _overlap1D(exp1, z1, n1, exp2, z2, n2);
    REAL ret = prefactor * sx * sy * sz;
    return ret;
}

REAL
overlap_PrimitiveGauss(
        const REAL exp1, const int l1, const int m1, const int n1, 
        const REAL x1,  const REAL y1, const REAL z1, const REAL norm_factor1,
        const REAL exp2, const int l2, const int m2, const int n2, 
        const REAL x2,  const REAL y2, const REAL z2, const REAL norm_factor2)
{
    REAL s = 0;
    s = _overlap3D( 
            exp1, l1, m1, n1, x1, y1, z1,
            exp2, l2, m2, n2, x2, y2, z2 );
    s *= norm_factor1 * norm_factor2;
    return s;
}

REAL
kinetic_PrimitiveGauss(
        const REAL exp1, const int l1, const int m1, const int n1, 
        const REAL x1,  const REAL y1, const REAL z1, const REAL norm_factor1,
        const REAL exp2, const int l2, const int m2, const int n2, 
        const REAL x2,  const REAL y2, const REAL z2, const REAL norm_factor2)
{
    REAL term1 = exp2 * (2*(l2+m2+n2)+3) * 
        _overlap3D(exp1,l1,m1,n1,x1,y1,z1,
                   exp2,l2,m2,n2,x2,y2,z2);

    REAL term2 = -2.0 * std::pow(exp2, 2) * (
            _overlap3D(exp1, l1, m1, n1, x1, y1, z1, exp2, l2 + 2, m2, n2, x2, y2, z2) +
            _overlap3D(exp1, l1, m1, n1, x1, y1, z1, exp2, l2, m2 + 2, n2, x2, y2, z2) + 
            _overlap3D(exp1, l1, m1, n1, x1, y1, z1, exp2, l2, m2, n2 + 2, x2, y2, z2) );

    REAL term3 = -0.5 * (
            l2*(l2-1) * _overlap3D(exp1, l1,m1,n1, x1,y1,z1, exp2, l2-2,m2,n2, x2,y2,z2) + 
            m2*(m2-1) * _overlap3D(exp1, l1,m1,n1, x1,y1,z1, exp2, l2,m2-2,n2, x2,y2,z2) + 
            n2*(n2-1) * _overlap3D(exp1, l1,m1,n1, x1,y1,z1, exp2, l2,m2,n2-2, x2,y2,z2) );

    return (term1 + term2 + term3) * norm_factor1 * norm_factor2;
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

REAL nuclear_attraction_PrimitiveGauss(
        const REAL exp1, const int l1, const int m1, const int n1, 
        const REAL x1, const REAL y1, const REAL z1, const REAL norm_factor1,
        const REAL exp2, const int l2, const int m2, const int n2, 
        const REAL x2, const REAL y2, const REAL z2, const REAL norm_factor2,
        const REAL xNuc, const REAL yNuc, const REAL zNuc)
{
    // <lhs | 1/Rzx | rhs>
    REAL gamma = exp1 + exp2;
    REAL dist2 = squaredNorm(x1, y1, z1, x2, y2, z2);

    REAL Rp_x = (exp1 * x1 + exp2 * x2)/gamma;
    REAL Rp_y = (exp1 * y1 + exp2 * y2)/gamma;
    REAL Rp_z = (exp1 * z1 + exp2 * z2)/gamma;

    REAL pc2 = squaredNorm(xNuc, yNuc, zNuc, Rp_x, Rp_y, Rp_z);
    REAL pre = 2.0 * M_PI / gamma * std::exp(-exp1 * exp2 * dist2 / gamma);

    g_list_array a_x, a_y, a_z;
    _g_list(l1, l2, Rp_x - x1, Rp_x - x2, Rp_x - xNuc, gamma, a_x);
    _g_list(m1, m2, Rp_y - y1, Rp_y - y2, Rp_y - yNuc, gamma, a_y);
    _g_list(n1, n2, Rp_z - z1, Rp_z - z2, Rp_z - zNuc, gamma, a_z);
    REAL s = 0.;
    for(int I = 0; I < 1+l1+l2; I++) {
        for(int J = 0; J < 1+m1+m2; J++) {
            for(int K = 0; K < 1+n1+n2; K++) {
                s += a_x[I] * a_y[J] * a_z[K] * boys(I+J+K, gamma*pc2);
            }
        }
    }
    s = s * pre * norm_factor1 * norm_factor2;
    return s;

}


typedef boost::array<REAL, 4*MAX_L + 1> c_list_array;

// Returns the C[I] in the paper pp.2320.
static void
_c_list(REAL exp1, REAL x1, int l1, REAL exp2, REAL x2, int l2, REAL px,
    REAL exp3, REAL x3, int l3, REAL exp4, REAL x4, int l4, REAL qx, c_list_array &c_list_out)
{
    std::fill(c_list_out.begin(), c_list_out.end(), 0.);
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


REAL electron_repulsion_PrimitiveGauss(
        const REAL exp1, const int l1, const int m1, const int n1, 
        const REAL x1, const REAL y1, const REAL z1, const REAL norm_factor1,
        const REAL exp2, const int l2, const int m2, const int n2, 
        const REAL x2, const REAL y2, const REAL z2, const REAL norm_factor2, 
        const REAL exp3, const int l3, const int m3, const int n3, 
        const REAL x3, const REAL y3, const REAL z3, const REAL norm_factor3,
        const REAL exp4, const int l4, const int m4, const int n4, 
        const REAL x4, const REAL y4, const REAL z4, const REAL norm_factor4 )
{
    REAL gamma1 = exp1 + exp2;
    REAL gamma2 = exp3 + exp4;

    REAL p_x = (exp1*x1 + exp2*x2)/gamma1;
    REAL p_y = (exp1*y1 + exp2*y2)/gamma1;
    REAL p_z = (exp1*z1 + exp2*z2)/gamma1;
    REAL q_x = (exp3*x3 + exp4*x4)/gamma2;
    REAL q_y = (exp3*y3 + exp4*y4)/gamma2;
    REAL q_z = (exp3*z3 + exp4*z4)/gamma2;

    REAL rab2 = squaredNorm(x1,y1,z1,x2,y2,z2);
    REAL rcd2 = squaredNorm(x3,y3,z3,x4,y4,z4);
    REAL rpq2 = squaredNorm(p_x, p_y, p_z, q_x, q_y, q_z);

    REAL delta = (1./gamma1 + 1./gamma2) / 4.;
    REAL s(0.);
    c_list_array c_list_x, c_list_y, c_list_z;
    _c_list(exp1, x1, l1, exp2, x2, l2, p_x, 
            exp3, x3, l3, exp4, x4, l4, q_x, c_list_x );
    _c_list(exp1, y1, m1, exp2, y2, m2, p_y, 
            exp3, y3, m3, exp4, y4, m4, q_y, c_list_y );
    _c_list(exp1, z1, n1, exp2, z2, n2, p_z, 
            exp3, z3, n3, exp4, z4, n4, q_z, c_list_z );

    for(int I = 0; I < 1 + l1 + l2 + l3 + l4; I++) {
        for(int J = 0; J < 1 + m1 + m2 + m3 + m4; J++) {
            for(int K = 0; K < 1 + n1 + n2 + n3 + n4; K++) {
                s += c_list_x[I] * c_list_y[J] * c_list_z[K] * boys(I+J+K, rpq2/(4.*delta));
            }
        }
    }
    s *= 2.0*std::pow(M_PI, 2.5)/ (gamma1 * gamma2)/ std::pow(gamma1 + gamma2, 0.5) * \
        std::exp(-(exp1*exp2*rab2/gamma1) - (exp3*exp4*rcd2/gamma2));

    s *= norm_factor1 * norm_factor2 * norm_factor3 * norm_factor4;
    return s;
}

}
