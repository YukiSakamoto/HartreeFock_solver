#define BOOST_TEST_MODULE "Overlap_Test"
#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include "../common.hpp"
#include "../gto.hpp"
#include "../gto_eval.hpp"
#include "test_common.hpp"


BOOST_AUTO_TEST_CASE(Hydrogen_R14_sto3g_integrals)
{
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

    Atom H1(1, 0., 0., 0.);
    Atom H2(1, 1.4, 0., 0.);

    // Overlap Integrals
    BOOST_CHECK(numerical_check( overlap_CGTO(h1, h1), 1.0, 0.0001 ));
    BOOST_CHECK(numerical_check( overlap_CGTO(h2, h2), 1.0, 0.0001 ));
    BOOST_CHECK(numerical_check( overlap_CGTO(h1, h2), 0.659318, 0.0001 ));
    BOOST_CHECK(numerical_check( overlap_CGTO(h2, h1), 0.659318, 0.0001 ));

    // Kinetic Integrals
    BOOST_CHECK(numerical_check( kinetic_CGTO(h1, h1), 0.7600, 0.0001 ));
    BOOST_CHECK(numerical_check( kinetic_CGTO(h2, h2), 0.7600, 0.0001 ));
    BOOST_CHECK(numerical_check( kinetic_CGTO(h1, h2), 0.2365, 0.0001 ));
    BOOST_CHECK(numerical_check( kinetic_CGTO(h2, h1), 0.2365, 0.0001 ));

    // Nuclear Attraction Integral for V1.
    BOOST_CHECK(numerical_check( nuclear_attraction_CGTO(h1, h1, H1), -1.2266, 0.0001 ));
    BOOST_CHECK(numerical_check( nuclear_attraction_CGTO(h2, h2, H2), -1.2266, 0.0001 ));

    BOOST_CHECK(numerical_check( electron_repulsion_CGTO(h1, h1, h1, h1), 0.7746, 0.0001 ) );
    BOOST_CHECK(numerical_check( electron_repulsion_CGTO(h2, h2, h2, h2), 0.7746, 0.0001 ) );
    BOOST_CHECK(numerical_check( electron_repulsion_CGTO(h1, h1, h2, h2), 0.5697, 0.0001 ) );
    BOOST_CHECK(numerical_check( electron_repulsion_CGTO(h2, h1, h1, h1), 0.4441, 0.0001 ) );
    BOOST_CHECK(numerical_check( electron_repulsion_CGTO(h2, h2, h2, h1), 0.4441, 0.0001 ) );
    BOOST_CHECK(numerical_check( electron_repulsion_CGTO(h2, h1, h2, h1), 0.2970, 0.0001 ) );
}

BOOST_AUTO_TEST_CASE(SP_same_center_overlap)
{
    ContractedGTO h1(0, 0, 0, 0., 0., 0.);
    h1.add_primitiveGTO(0.444635, 0.168856);
    h1.add_primitiveGTO(0.535328, 0.623913);
    h1.add_primitiveGTO(0.154329, 3.42525);
    h1.normalize();

    ContractedGTO h_2px(1, 0, 0, 0., 0., 0.);
    h_2px.add_primitiveGTO(0.444635, 0.168856);
    h_2px.add_primitiveGTO(0.535328, 0.623913);
    h_2px.add_primitiveGTO(0.154329, 3.42525);
    h_2px.normalize();

    ContractedGTO h_2py(0, 1, 0, 0., 0., 0.);
    h_2py.add_primitiveGTO(0.444635, 0.168856);
    h_2py.add_primitiveGTO(0.535328, 0.623913);
    h_2py.add_primitiveGTO(0.154329, 3.42525);
    h_2py.normalize();

    ContractedGTO h_2pz(0, 0, 1, 0., 0., 0.);
    h_2pz.add_primitiveGTO(0.444635, 0.168856);
    h_2pz.add_primitiveGTO(0.535328, 0.623913);
    h_2pz.add_primitiveGTO(0.154329, 3.42525);
    h_2pz.normalize();

    BOOST_CHECK( numerical_check(overlap_CGTO(h1, h_2px), 0., 0.0001));
    BOOST_CHECK( numerical_check(overlap_CGTO(h1, h_2py), 0., 0.0001));
    BOOST_CHECK( numerical_check(overlap_CGTO(h1, h_2pz), 0., 0.0001));
}

BOOST_AUTO_TEST_CASE(SP_different_center_overlap)
{
    ContractedGTO h1(0, 0, 0, 0., 0., 0.);
    h1.add_primitiveGTO(0.444635, 0.168856);
    h1.add_primitiveGTO(0.535328, 0.623913);
    h1.add_primitiveGTO(0.154329, 3.42525);
    h1.normalize();

    ContractedGTO h_2px(1, 0, 0, 1.4, 0., 0.);
    h_2px.add_primitiveGTO(0.444635, 0.168856);
    h_2px.add_primitiveGTO(0.535328, 0.623913);
    h_2px.add_primitiveGTO(0.154329, 3.42525);
    h_2px.normalize();

    ContractedGTO h_2py(0, 1, 0, 1.4, 0., 0.);
    h_2py.add_primitiveGTO(0.444635, 0.168856);
    h_2py.add_primitiveGTO(0.535328, 0.623913);
    h_2py.add_primitiveGTO(0.154329, 3.42525);
    h_2py.normalize();

    ContractedGTO h_2pz(0, 0, 1, 1.4, 0., 0.);
    h_2pz.add_primitiveGTO(0.444635, 0.168856);
    h_2pz.add_primitiveGTO(0.535328, 0.623913);
    h_2pz.add_primitiveGTO(0.154329, 3.42525);
    h_2pz.normalize();

    //BOOST_CHECK( numerical_check(overlap_CGTO(h1, h_2px), 0., 0.0001));
    std::cout << "<S|Px>: " << overlap_CGTO(h1, h_2px) << std::endl;
    BOOST_CHECK( overlap_CGTO(h1, h_2px) < 0. );
    BOOST_CHECK( numerical_check(overlap_CGTO(h1, h_2py), 0., 0.0001));
    BOOST_CHECK( numerical_check(overlap_CGTO(h1, h_2pz), 0., 0.0001));
}


