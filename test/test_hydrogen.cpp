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

using namespace MOSolver;

BOOST_AUTO_TEST_CASE(Hydrogen_1S_system)
{
    Atom H1(1, 0., 0., 0.);
    Atom H2(1, 1.4, 0., 0.);
    System H2_molecule;
    H2_molecule.add_atom(H1);
    H2_molecule.add_atom(H2);

    BOOST_CHECK( numerical_check(H2_molecule.total_electrons(), 2, 0.00) );
    BOOST_CHECK( numerical_check(H2_molecule.nspin_, 0, 0.00) );

    ContractedGTO h1(0, 0, 0, 0., 0., 0.);
    h1.add_primitiveGTO(0.444635, 0.168856);
    h1.add_primitiveGTO(0.535328, 0.623913);
    h1.add_primitiveGTO(0.154329, 3.42525);
    h1.normalize();

    ContractedGTO h2( 0, 0, 0, 1.4, 0., 0.);
    h2.add_primitiveGTO(0.444635, 0.168856);
    h2.add_primitiveGTO(0.535328, 0.623913);
    h2.add_primitiveGTO(0.154329, 3.42525);
    h2.normalize();

    CGTOs bfs;
    bfs.push_back(h1);
    bfs.push_back(h2);

    REAL NEI = H2_molecule.nuclear_repulsion();
    //BOOST_CHECK( numerical_check(NEI, 1.8310-1.1167, 0.001) );
    BOOST_CHECK( numerical_check(NEI, 0.7143, 0.001) );
    BOOST_CHECK( std::abs(rhf(bfs, H2_molecule) - uhf(bfs, H2_molecule)) < 0.001 );
}

