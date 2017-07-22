
#define _USE_MATH_DEFINES
#include <cmath>
#include "Eigen/Core"
#include "common.hpp"
#include "system.hpp"
#include "orbitals.hpp"
#include "grid.hpp"
#include <cstdio>
#include <string>
#include <vector>


#include "lebedev_194.h"

struct radial_grid legendre32[] = {
	{ -0.997263861849 , 0.007018610009  }, 
	{ -0.985611511545 , 0.016274394731  }, 
	{ -0.964762255588 , 0.025392065309  }, 
	{ -0.934906075938 , 0.034273862913  }, 
	{ -0.896321155766 , 0.042835898022  }, 
	{ -0.849367613733 , 0.050998059262  }, 
	{ -0.794483795968 , 0.058684093479  }, 
	{ -0.732182118740 , 0.065822222776  }, 
	{ -0.663044266930 , 0.072345794109  }, 
	{ -0.587715757241 , 0.078193895787  }, 
	{ -0.506899908932 , 0.083311924227  }, 
	{ -0.421351276131 , 0.087652093004  }, 
	{ -0.331868602282 , 0.091173878696  }, 
	{ -0.239287362252 , 0.093844399081  }, 
	{ -0.144471961583 , 0.095638720079  }, 
	{ -0.048307665688 , 0.096540088515  }, 
	{  0.048307665688 , 0.096540088515  }, 
	{  0.144471961583 , 0.095638720079  }, 
	{  0.239287362252 , 0.093844399081  }, 
	{  0.331868602282 , 0.091173878696  }, 
	{  0.421351276131 , 0.087652093004  }, 
	{  0.506899908932 , 0.083311924227  }, 
	{  0.587715757241 , 0.078193895787  }, 
	{  0.663044266930 , 0.072345794109  }, 
	{  0.732182118740 , 0.065822222776  }, 
	{  0.794483795968 , 0.058684093479  }, 
	{  0.849367613733 , 0.050998059262  }, 
	{  0.896321155766 , 0.042835898022  }, 
	{  0.934906075938 , 0.034273862913  }, 
	{  0.964762255588 , 0.025392065309  }, 
	{  0.985611511545 , 0.016274394731  }, 
	{  0.997263861849 , 0.007018610009  }, 
};

// Referencing Papers:
//   A.D.Becke,  J.Chem.Phys 1988, 88, 2547-2553
//************************************************************
//
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

REAL weight(const Vector3Real &pos, const Vector3Real &atom_pos, const Vector3Real &another_atom_pos)
{
    REAL Rij = (atom_pos - another_atom_pos).norm();
    REAL lambda1   = (pos- atom_pos).norm();
    REAL lambda2   = (pos- another_atom_pos).norm();
    REAL u = (lambda1 - lambda2)/Rij;
    return s_n(3, u);
}

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
// Numerical Quadrature Points
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

VectorXReal
LegendreRoots(int n)
{
    if (n != 32) {  throw;  }
    VectorXReal ret = VectorXReal::Zero(n);
    for(int i = 0; i < n; i++) {
        ret[i] = legendre32[i].x;
    }
    return ret;
}
VectorXReal
LegendreWeights(int n)
{
    if (n != 32) {  throw;  }
    VectorXReal ret = VectorXReal::Zero(n);
    for(int i = 0; i < n; i++) {
        ret[i] = legendre32[i].weight;
    }
    return ret;

}

//************************************************************
// Grid Related 
//************************************************************
void
PrintGridWeight(const AtomicGrid &grid)
{
    int ngrid = grid.size();
    for(int i = 0; i < ngrid; i++) {
        std::printf("  %18.10f\t%18.10f\t%18.10f\t%18.10f\t%18.10f\t%18.10f\tr = %18.10f\n", 
                grid.grid_points[i].pos[0], grid.grid_points[i].pos[1], grid.grid_points[i].pos[2],
                grid.grid_points[i].weight1, grid.grid_points[i].weight2,
                grid.grid_points[i].weight1 * grid.grid_points[i].weight2, grid.grid_points[i].radial_);
    }
}
//
//
REAL
calc_BeckeWeight(const GridPoint &grid_point, const System& system, const int index)
{
    int natom = system.size();
    REAL w_this = 0.;
    REAL w_total = 0.;
    for(int i = 0; i < natom; i++) {
        REAL w = 1.;
        for(int j = 0; j < natom; j++) {
            if (i == j) { continue; }
            w *= weight(grid_point.pos, system[i].center, system[j].center);    
        }
        if (i == index) {   w_this = w;    }
        w_total += w;
    }
    return w_this/w_total;
}

REAL
getRmax(int atomic_number)
{
    // XXX Note:
    //  this function returns the Rmax in Bohr Unit, NOT in Angstrom.
    //
    // XXX For now, C and H ONLY.
    REAL r = 0.;
    if (atomic_number == 1) {
        r = 0.25;
    } else if (atomic_number == 6) {
        r = 0.70 / 2.;
    } else {
        ;
    }
    return Angstrom2Bohr(r);
}

AtomicGrid
GenerateAtomicGrid(const System &system, const int index, int nrad = 32, int nang = 194) 
{
    AtomicGrid atomic_grid(system[index].center);
    REAL Rmax = getRmax(system[index].atomic_number);

    std::printf("Generating Grid of Atom No. %03d\n", index);
    // 1.  determine the radial quadrature point
    VectorXReal x_quadpoints = LegendreRoots(nrad);
    VectorXReal x_weights    = LegendreWeights(nrad);

    // 2. map the quadrature points to the real space.
    //    each elements corresponds to the r_i.
    VectorXReal r_quadpoints = coordinate_map(x_quadpoints, Rmax);
    
    std::cout << "radial grid determined" << std::endl;

    for(int i = 0; i < nrad; i++) {
        REAL x = x_quadpoints[i];
        REAL r = r_quadpoints[i];
        // Generate Spherical Grid
        for(int j = 0; j < nang; j++) {
            Vector3Real lebedev_point = Vector3Real(lebedev_grid_194[j].x, lebedev_grid_194[j].y, lebedev_grid_194[j].z);
            REAL lebedev_weight = lebedev_grid_194[j].weight;

            REAL dr = Rmax*2/std::pow(1-x,2);   // <--- dr/dx 
            REAL weight = 4*M_PI*r*r*dr*lebedev_weight;
            GridPoint gp(system[index].center + lebedev_point*r, weight);
            gp.weight2 = calc_BeckeWeight(gp, system, index);
            gp.radial_ = r;
            atomic_grid.append(gp);
        }
    }
    return atomic_grid;
}

std::vector<AtomicGrid>
GenerateGrid(const System &system)
{
    std::vector<AtomicGrid> grid;
    int natoms = system.size();
    std::printf("natoms: %d\n", natoms);
    size_t ngrid = 0;
    for(int i = 0; i < natoms; i++) {
        grid.push_back( GenerateAtomicGrid(system, i) );
        ngrid += grid[i].size();
    }
    std::cout << ngrid << " grids generated" << std::endl;
    return grid;
}

