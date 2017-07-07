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

    CGTOs bfs;
    {   // Orbital Generation
        std::vector<REAL> C_1s_alpha;
        std::vector<REAL> C_1s_coeff;
        {
            C_1s_alpha.push_back( 71.6168370 );
            C_1s_alpha.push_back( 13.0450960 );
            C_1s_alpha.push_back(  3.5305122 );
            C_1s_coeff.push_back( 0.15432897  );
            C_1s_coeff.push_back( 0.53532814  );
            C_1s_coeff.push_back( 0.44463454  );
        }
        bfs.add_orbitals(0, C.center, C_1s_alpha, C_1s_coeff);

        std::vector<REAL> C_2sp_alpha;
        std::vector<REAL> C_2s_coeff;
        std::vector<REAL> C_2p_coeff;
        {
            C_2sp_alpha.push_back(  2.9412494 );
            C_2sp_alpha.push_back(  0.6834831 );
            C_2sp_alpha.push_back(  0.2222899 );
            C_2s_coeff.push_back(-0.09996723 );
            C_2s_coeff.push_back( 0.39951283 );
            C_2s_coeff.push_back( 0.70011547 );
            C_2p_coeff.push_back( 0.15591627 );
            C_2p_coeff.push_back( 0.60768372 );
            C_2p_coeff.push_back( 0.39195739 );
        }
        bfs.add_orbitals(0, C.center, C_2sp_alpha, C_2s_coeff);
        bfs.add_orbitals(1, C.center, C_2sp_alpha, C_2p_coeff);

        std::vector<REAL> H_1s_alpha;
        std::vector<REAL> H_1s_coeff;
        {
            H_1s_coeff.push_back(0.444635);
            H_1s_coeff.push_back(0.535328);
            H_1s_coeff.push_back(0.154329);
            H_1s_alpha.push_back(0.168856);
            H_1s_alpha.push_back(0.623913);
            H_1s_alpha.push_back(3.42525 );
        }
        bfs.add_orbitals(0, H1.center, H_1s_alpha, H_1s_coeff);
        bfs.add_orbitals(0, H2.center, H_1s_alpha, H_1s_coeff);
        bfs.add_orbitals(0, H3.center, H_1s_alpha, H_1s_coeff);
        bfs.add_orbitals(0, H4.center, H_1s_alpha, H_1s_coeff);
    }

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
