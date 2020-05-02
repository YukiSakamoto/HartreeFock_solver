#pragma once

namespace libGaussEval {

typedef double REAL;
const int MAX_L = 2;    

REAL overlap_PrimitiveGauss(
        const REAL exp1, const int l1, const int m1, const int n1, 
        const REAL x1, const REAL y1, const REAL z1, const REAL norm_factor1,
        const REAL exp2, const int l2, const int m2, const int n2, 
        const REAL x2, const REAL y2, const REAL z2, const REAL norm_factor2 );

REAL kinetic_PrimitiveGauss(
        const REAL exp1, const int l1, const int m1, const int n1, 
        const REAL x1, const REAL y1, const REAL z1, const REAL norm_factor1,
        const REAL exp2, const int l2, const int m2, const int n2, 
        const REAL x2, const REAL y2, const REAL z2, const REAL norm_factor2 );

REAL nuclear_attraction_PrimitiveGauss(
        const REAL exp1, const int l1, const int m1, const int n1, 
        const REAL x1, const REAL y1, const REAL z1, const REAL norm_factor1,
        const REAL exp2, const int l2, const int m2, const int n2, 
        const REAL x2, const REAL y2, const REAL z2, const REAL norm_factor2,
        const REAL xNuc, const REAL yNuc, const REAL zNuc );

REAL electron_repulsion_PrimitiveGauss(
        const REAL exp1, const int l1, const int m1, const int n1, 
        const REAL x1, const REAL y1, const REAL z1, const REAL norm_factor1,
        const REAL exp2, const int l2, const int m2, const int n2, 
        const REAL x2, const REAL y2, const REAL z2, const REAL norm_factor2, 
        const REAL exp3, const int l3, const int m3, const int n3, 
        const REAL x3, const REAL y3, const REAL z3, const REAL norm_factor3,
        const REAL exp4, const int l4, const int m4, const int n4, 
        const REAL x4, const REAL y4, const REAL z4, const REAL norm_factor4 );


} // namespace LibGaussEval
