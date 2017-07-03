
#include "gto.hpp"
#include <cmath>

REAL 
PrimitiveGTO::value(REAL rx, REAL ry, REAL rz) {
    Vector3Real r(Vector3Real(rx, ry, rz));
    return value(r);
}

REAL
PrimitiveGTO::value(Vector3Real v) {
    Vector3Real d(v - center);
    REAL norm2 = d.squaredNorm();
    REAL phase = std::pow(v[0] - center[0], l) * std::pow(v[1] - center[1], m) * std::pow(v[2] - center[2], n);
    return phase * std::exp(-1 * exponent * norm2) * norm_factor;
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

