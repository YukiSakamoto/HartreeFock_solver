#define BOOST_TEST_MODULE "H2_test"
#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include "../common.hpp"
#include "../gto.hpp"
#include "../gto_eval.hpp"
#include "../hf.hpp"
#include "test_common.hpp"



BOOST_AUTO_TEST_CASE(CH4_rhf_system)
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

    BOOST_CHECK( numerical_check(CH4.total_electrons(), 10, 0.00) );
    BOOST_CHECK( numerical_check(CH4.nspin() , 0, 0.00) );

    CGTOs bfs = generate_bfs(CH4, "sto3g.dat");
    REAL rhf_en = rhf(bfs, CH4);
    std::cerr << rhf_en << std::endl;
    REAL uhf_en = uhf(bfs, CH4);
    std::cerr << uhf_en << std::endl;
    const REAL Etot = -39.7269;

    //BOOST_CHECK( numerical_check(NEI, 1.8310-1.1167, 0.001) );
    BOOST_CHECK( std::abs(rhf_en - Etot) < 0.001 );
    BOOST_CHECK( std::abs(rhf_en - uhf_en) < 0.001 );
    
}

