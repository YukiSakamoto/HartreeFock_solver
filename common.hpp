#pragma once
#include "Eigen/Core"

// numerical precision: Double 
typedef double REAL;
typedef Eigen::VectorXd VectorXReal;
typedef Eigen::Vector3d Vector3Real;
typedef Eigen::Vector2d Vector2Real;
typedef Eigen::MatrixXd MatrixXReal;

// Constant Values
static const REAL factor_Bohr2Angstrom = 0.529177249;
static const REAL factor_Angstrom2Bohr = 1.8897259885789;

inline REAL
Bohr2Angstrom(REAL x)
{
    return x*factor_Bohr2Angstrom;
}

inline Vector3Real
Bohr2Angstrom(Vector3Real r)
{
    return r * factor_Bohr2Angstrom;
}


inline REAL
Angstrom2Bohr(REAL x)
{
    return x * factor_Angstrom2Bohr;
}

inline Vector3Real
Angstrom2Bohr(Vector3Real r)
{
    return r * factor_Angstrom2Bohr;
}

inline Vector3Real
Angstrom2Bohr(REAL x, REAL y, REAL z)
{
    Vector3Real r = Vector3Real(x,y,z);
    return r * factor_Angstrom2Bohr;
}
