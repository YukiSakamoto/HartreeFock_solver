
#include "gto.hpp"
#include <cmath>

REAL 
PrimitiveGTO::value(REAL rx, REAL ry, REAL rz) {
    REAL dx = rx - x;
    REAL dy = ry - y;
    REAL dz = rz - z;
    REAL norm2 = dx*dx + dy*dy + dz*dz;
    return std::exp(-1 * exponent * norm2) * norm_factor;
}

REAL
PrimitiveGTO::value(Vector3Real v) {
    Vector3Real d(v - center);
    REAL norm2 = d.squaredNorm();
    return std::exp(-1 * exponent * norm2) * norm_factor;
}

void
PrimitiveGTO::normalize() {
    REAL numerator = std::pow(2, 2*(l+m+n)+1.5) * std::pow(exponent, l+m+n+1.5);
    REAL denominator = factorial2(2*l-1) * factorial2(2*m-1) * factorial2(2*n-1) * std::pow(M_PI, 1.5);

    if (denominator == REAL(0.)) {
        std::cout << __LINE__ << std::endl;
        throw;
    }
    norm_factor = std::sqrt(numerator / denominator);
}

