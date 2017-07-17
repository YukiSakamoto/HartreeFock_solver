
#define _USE_MATH_DEFINES
#include <cmath>
#include "Eigen/Core"
#include "common.hpp"
#include <cstdio>
#include <string>
#include <vector>

#include "lebedev_194.h"

// Referencing Papers:
//   A.D.Becke,  J.Chem.Phys 1988, 88, 2547-2553

//************************************************************
// Becke's Functions
// Becke 1988 (eq. 19)
inline 
REAL p(REAL u)
{
    return (3./2.) * u - (1./2.) * std::pow(u,3);
}

// Becke 1988 (eq. 20)
inline
REAL f_n(size_t n, REAL u)
{
    if (n <= 0) {   throw;  }
    REAL val = u;
    for(int i = 0; i < n; i++) {    
        val = p(val); 
    }
    return val;
}
// Becke 1988 (eq. 21)
REAL s_n(size_t n, REAL u)
{   return (1./2.) * (1-f_n(n,u));    }

// Becke 1988 (eq. 25)
REAL coordinate_map(REAL x, const REAL Rmax)
{
    // Map the integration interval 
    //   -1 <= x < 1 to 0<= r < inf to 
    return Rmax * (1.+x)/(1.-x);

    // XXX Note:  dr/dx = Rmax*2x/(1-x)^2
}

VectorXReal
coordinate_map(const VectorXReal &X, const REAL Rmax) 
{
    int len = X.size();
    VectorXReal ret = VectorXReal::Zero(len);
    for(int i = 0; i < len; i++) {
        ret[i] = coordinate_map(X[i], Rmax);
    }
    return ret;
}
//************************************************************
// Numerical Quadrature 
//************************************************************
VectorXReal
Chebyshev2ndRoots(int n)
{
    if (n < 0) {    throw;  }
    VectorXReal ret = VectorXReal::Zero(n);
    for(int i = 0; i < n; i++) {
        // abramowitz and stegun 25.4.40
        ret[i] = std::cos(M_PI*(i+1)/(n+1));
    }
    return ret;
}

VectorXReal
Chebyshev2ndWeights(int n)
{
    if (n < 0) {    throw;  }
    VectorXReal ret = VectorXReal::Zero(n);
    for(int i = 0; i < n; i++) {
        ret[i] = M_PI/(n+1)*std::pow( std::sin(M_PI*(i+1)/(n+1)), 2 );
    }
    return ret;
}

//************************************************************
// Grid Related 
//************************************************************
struct GridPoint {
    REAL radial_weight;
    REAL spherical_weight;
    Vector3Real pos;
    GridPoint(const Vector3Real &x, REAL weight = 0.): 
        pos(x), radial_weight(weight)
    {}
};

struct AtomicGrid {
    std::vector<GridPoint> grid_points;
    Vector3Real pos;
    
    void push_back(const GridPoint& gp) {
        this->grid_points.push_back(gp);
    }
    REAL integrate() const 
    {
        REAL sum = 0.;
        size_t length = this->grid_points.size(); 
        for(size_t i = 0; i < length; i++) {
            const GridPoint &grid = grid_points[i];
        }
        return sum;
    }
};

//AtomicGrid
void
GenerateAtomicGrid(const Vector3Real &atom_pos, REAL Rmax, int nrad = 70, int nang = 194) 
{
    // 1.  determine the radial quadrature point
    VectorXReal x_quadpoints = Chebyshev2ndRoots(nrad);
    VectorXReal x_weights    = Chebyshev2ndWeights(nrad);

    // 2. map the quadrature points to the real space.
    //    each elements corresponds to the r_i.
    VectorXReal r_quadpoints = coordinate_map(x_quadpoints, Rmax);

    for(int i = 0; i < nrad; i++) {
        REAL r = r_quadpoints[i];
        // Generate Spherical Grid
        Vector3Real lebedev_point = Vector3Real(lebedev_grid_194[i].x, lebedev_grid_194[i].y, lebedev_grid_194[i].z);
        REAL lebedev_weight = lebedev_grid_194[i].weight;
        REAL weight = 4*M_PI*r*r * lebedev_weight;
        GridPoint(atom_pos + lebedev_point * r);
    }

    // memorize
}


REAL f(REAL x, REAL y, REAL z) 
{
    // Vexact = 
    return 1. + x + std::pow(y,2) + std::pow(x,2)*y + std::pow(x,4) + std::pow(y,5) + std::pow(x*y*z, 2);
}

REAL SphereVolume(REAL R)
{
    REAL s = 0.;
    REAL Sexact = 4.*M_PI*R*R*R/3.;
    REAL dr = 0.00001;
    for(REAL radius = 0.; radius <= R; radius += dr) {
        REAL S_r = 4*M_PI*radius*radius;
        for(int i = 0; i < sizeof(lebedev_grid_194)/sizeof(lebedev_grid); i++) {
            REAL weight = lebedev_grid_194[i].weight;
            s += S_r * weight * dr;
        }
    }
    std::printf("Quadrature: %.8f, Vexact: %.8f\n", s, Sexact);
    return s;
}

REAL g(REAL x)
{
    return 4./(10+x*x);
}

int main(void)
{
    // integrate 0->1: r^2
    //REAL s = 0.;
    //REAL S = 4*M_PI;
    //for(int i = 0; i < sizeof(lebedev_grid_194)/sizeof(lebedev_grid); i++) {
    //    Vector3Real r = Vector3Real(lebedev_grid_194[i].x, lebedev_grid_194[i].y, lebedev_grid_194[i].z);
    //    REAL weight = lebedev_grid_194[i].weight;
    //    s += f(r[0], r[1], r[2]) * weight * S;
    //}
    //REAL Vexact = 216*M_PI/35.;
    //std::printf("Vexact: %.8f,  Vquadrature:  %.8f\n", Vexact, s);
    //
    int nrad = 40;
    REAL s = 0;
    VectorXReal xquadpoints = Chebyshev2ndRoots(nrad);
    VectorXReal xweights    = Chebyshev2ndWeights(nrad);
    
    for(int i = 0; i < xquadpoints.size(); i++) {
        std::printf("x%02d: %12.8f\tw%02d: %12.8f\n", i, xquadpoints[i], i, xweights[i]);
        s += g(xquadpoints[i]) * xweights[i];
    }
    std::printf("Chevyshev Quadrature of g(x): %.8f\n", s);
    //SphereVolume(2.0);
    int n = 10000;
    REAL dx = 2./n;
    s = 0.;
    for(REAL x = -1.; x <= 1.;  x += dx) {
        s += g(x)*dx;
    }
    std::printf("Sum of g(x)dx (dx = %.8f) : %.8f\n", dx, s);
    return 0;
}
