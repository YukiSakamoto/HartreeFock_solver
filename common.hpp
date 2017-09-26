#pragma once
#include "Eigen/Core"

// numerical precision: Double 
typedef double REAL;
typedef Eigen::VectorXd VectorXReal;
typedef Eigen::Vector3d Vector3Real;
typedef Eigen::Vector2d Vector2Real;
typedef Eigen::MatrixXd MatrixXReal;

// For the vector of magnetic quantum number
typedef Eigen::Vector3i Vector3Int;

// Maximum Angular Momentum
const int MAX_L = 2;    
// 
const std::string Elements[] = {
    "0", 
    "H", "He", 
    "Li", "Be", "B", "C", "N", "O", "F", "Ne", 
    "Na", "Mg", "Al", "Si", "P" , "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", 
    "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "",
};

inline int 
get_atomic_number(const std::string &element)
{
    for(int i = 0; Elements[i] != ""; i++) {
        if (Elements[i] == element) {
            return i;
        }
    }
    return -1;
}

inline
int angular_momentum(char type)
{
    char angular_list[] = {'S', 'P', 'D', 'F', 0};
    for(int i = 0; angular_list[i] != 0; i++) {
        if (angular_list[i] == type) {
            return i;
        }
    }
    return -1;
}


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
