#include "common.hpp"

#include "math_helper.hpp"
#include "gto.hpp"
#include "gto_eval.hpp"
#include "orbitals.hpp"
#include "system.hpp"
#include "hf.hpp"

#include <iostream>


REAL
CH4(REAL /*ignored*/)
{
    Atom C(6, Angstrom2Bohr(0., 0., 0.) );
    Atom H1(1,Angstrom2Bohr(  0.000000,    0.000000,    1.083010) );
    Atom H2(1,Angstrom2Bohr(  0.000000,    1.021071,   -0.361003) );
    Atom H3(1,Angstrom2Bohr(  0.884274,   -0.510536,   -0.361003) );
    Atom H4(1,Angstrom2Bohr( -0.884274,   -0.510536,   -0.361003) );
    System CH4;
    CH4.add_atom(C);
    CH4.add_atom(H1);
    CH4.add_atom(H2);
    CH4.add_atom(H3);
    CH4.add_atom(H4);

    CGTOs bfs = generate_bfs(CH4, "sto3g.dat");
    //CGTOs bfs = generate_bfs(CH4, "p631G.dat");
    rhf(bfs, CH4);

    return 0.;
}

int main(void)
{
    CH4(1.4);
    //HeH(1.4632);
    //He(1.4632);
    return 0;
}
