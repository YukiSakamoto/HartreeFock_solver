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
#include "../orbitals.hpp"
#include "../hf.hpp"
#include "test_common.hpp"

BOOST_AUTO_TEST_CASE(HeH_R14632_sto3g)
{
    REAL r12 = 1.4632;
    Atom H(1, 0. , 0., 0.);
    Atom He(2, r12, 0., 0.);
    System HeH;
    HeH.add_atom(H);
    HeH.add_atom(He);
    HeH.set_charge(+1);

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

    REAL Etot = rhf(bfs, HeH);
    BOOST_CHECK(numerical_check(Etot, -2.86066, 0.0001));
    BOOST_CHECK(numerical_check(uhf(bfs, HeH), -2.86066, 0.0001));
}
