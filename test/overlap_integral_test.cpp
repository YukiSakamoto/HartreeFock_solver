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


BOOST_AUTO_TEST_CASE(Overlap_H_1s)
{
    ContractedGTO h1(0, 0, 0, 0., 0., 0.);
    h1.add_primitiveGTO(0.444635, 0.168856);
    h1.add_primitiveGTO(0.535328, 0.623913);
    h1.add_primitiveGTO(0.154329, 3.42525);
    h1.normalize();

    ContractedGTO h3( 0, 0, 0, 1.4, 0., 0.);
    h3.add_primitiveGTO(0.444635, 0.168856);
    h3.add_primitiveGTO(0.535328, 0.623913);
    h3.add_primitiveGTO(0.154329, 3.42525);
    h3.normalize();

    BOOST_CHECK(numerical_check( overlap_CGTO(h1, h1), 1.0, 0.0001 ));
    BOOST_CHECK(numerical_check( overlap_CGTO(h3, h3), 1.0, 0.0001 ));
    BOOST_CHECK(numerical_check( overlap_CGTO(h1, h3), 0.659318, 0.0001 ));
    BOOST_CHECK(numerical_check( overlap_CGTO(h3, h1), 0.659318, 0.0001 ));
}

