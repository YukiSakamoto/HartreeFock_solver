#include "Eigen/Core"
#include "common.hpp"

#include "math_helper.hpp"
#include "gto.hpp"
#include "gto_eval.hpp"
#include "orbitals.hpp"
#include "system.hpp"
#include "hf.hpp"

#include <iostream>

#include <Eigen/Eigen> 

REAL
HeH(REAL r12)
{
    //==================================================
    //  Prepare the system
    //==================================================
    Atom H(1, 0. , 0., 0.);
    Atom He(2, r12, 0., 0.);
    System HeH;
    HeH.add_atom(H);
    HeH.add_atom(He);
    HeH.set_charge(+1);

    //==================================================
    //  Prepare the Basis Functions
    //==================================================
    CGTOs bfs;

    ContractedGTO he1s(0, 0, 0, r12, 0., 0.);
    he1s.add_primitiveGTO(0.444635, 0.480844);
    he1s.add_primitiveGTO(0.535328, 1.776691);
    he1s.add_primitiveGTO(0.154329, 9.753934);
    he1s.normalize();
    bfs.push_back(he1s);

    ContractedGTO h1s(0, 0, 0, 0., 0., 0.);
    h1s.add_primitiveGTO(0.444635, 0.168856);
    h1s.add_primitiveGTO(0.535328, 0.623913);
    h1s.add_primitiveGTO(0.154329, 3.42525);
    h1s.normalize();
    bfs.push_back(h1s);

    rhf(bfs, HeH);
    return 0.;
}

REAL 
H2(REAL r12)
{
    CGTOs bfs;

    ContractedGTO h1(0, 0, 0, 0., 0., 0.);
    h1.add_primitiveGTO(0.444635, 0.168856);
    h1.add_primitiveGTO(0.535328, 0.623913);
    h1.add_primitiveGTO(0.154329, 3.42525);
    h1.normalize();

    ContractedGTO h2( 0, 0, 0, r12, 0., 0.);
    h2.add_primitiveGTO(0.444635, 0.168856);
    h2.add_primitiveGTO(0.535328, 0.623913);
    h2.add_primitiveGTO(0.154329, 3.42525);
    h2.normalize();

    bfs.push_back(h1);
    bfs.push_back(h2);

    System H2_molecule;
    Atom H1(1, 0. , 0., 0.);
    Atom H2(1, r12, 0., 0.);
    H2_molecule.add_atom(H1);
    H2_molecule.add_atom(H2);

    rhf(bfs, H2_molecule);
    return 0.;
}

REAL
He(REAL /*ignored*/)
{
    Atom He_atom(2, 0, 0., 0.);
    System He_;
    He_.add_atom(He_atom);

    CGTOs bfs;

    ContractedGTO he1s(0, 0, 0, 0., 0., 0.);
    he1s.add_primitiveGTO(0.444635, 0.480844);
    he1s.add_primitiveGTO(0.535328, 1.776691);
    he1s.add_primitiveGTO(0.154329, 9.753934);
    he1s.normalize();
    bfs.push_back(he1s);
    rhf(bfs, He_);

    return 0.;
}


int main(void)
{
    //H2(1.4);
    HeH(1.4632);
    //He(1.4632);
    return 0;
}
