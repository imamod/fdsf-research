#include "Common.h"

#ifdef HIGH_PRECISION
TEST_CASE("BOOST_GAMMA") {
    using namespace boost::multiprecision;
    std::setprecision(std::numeric_limits<BmpReal>::max_digits10);
    // Operations at fixed precision and full numeric_limits support:
    cpp_dec_float_100 b = 2;
    std::cout << std::numeric_limits<cpp_dec_float_100>::digits << std::endl;
    // Note that digits10 is the same as digits, since we're base 10! :
    std::cout << std::numeric_limits<cpp_dec_float_100>::digits10 << std::endl;
    // We can use any C++ std lib function, lets print all the digits as well:
    std::cout << std::setprecision(std::numeric_limits<cpp_dec_float_100>::max_digits10)
        << log(b) << std::endl << log(BmpReal(2)) << std::endl; // print log(2)
                                                                 // We can also use any function from Boost.Math:
    std::cout << boost::math::tgamma(b) << std::endl;
    // These even work when the argument is an expression template:
    std::cout << boost::math::tgamma(b * b) << std::endl;
    // And since we have an extended exponent range we can generate some really large 
    // numbers here (4.0238726007709377354370243e+2564):
    std::cout << boost::math::tgamma(cpp_dec_float_100(1000)) << std::endl;
    std::cout << boost::math::tgamma(cpp_dec_float_100(3.0 / 2.0)) << std::endl;
}

TEST_CASE("Calc_PI") {
    BmpReal my_pi = boost::math::constants::pi<BmpReal>();
    std::cout << "my_pi = " << my_pi << std::endl;
    std::cout << "4atg1 = " << 4 * boost::multiprecision::atan(BmpReal(1)) << std::endl;
    std::cout << "I1(0) = " << my_pi*my_pi / 12 << std::endl;
}
#endif
