#define BOOST_TEST_MODULE "Overlap_Test"
#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include "../common.hpp"
#include "../math_helper.hpp"
#include "test_common.hpp"

using namespace MOSolver;

BOOST_AUTO_TEST_CASE(NetfreeModel_test_constructor)
{
    BOOST_CHECK( numerical_check( boys(0., 0.), 1.0, 0.001 ) );
}


