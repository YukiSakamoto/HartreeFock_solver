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

using namespace MOSolver;

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

    CGTOs bfs = generate_bfs(HeH, "p631G.dat");
    rhf(bfs, HeH);
    return 0.;
}

REAL 
H2(REAL r12)
{
    System H2_molecule;
    Atom H1(1, 0. , 0., 0.);
    Atom H2(1, r12, 0., 0.);
    H2_molecule.add_atom(H1);
    H2_molecule.add_atom(H2);

    CGTOs bfs = generate_bfs(H2_molecule, "sto3g.dat");
    rhf(bfs, H2_molecule);
    return 0.;
}

REAL
He(REAL /*ignored*/)
{
    Atom He_atom(2, 0, 0., 0.);
    System He_;
    He_.add_atom(He_atom);

    CGTOs bfs = generate_bfs(He_, "sto3g.dat");
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
