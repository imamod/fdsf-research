#define CATCH_CONFIG_MAIN
#include "TestCommon.h"

namespace {
    bmp_real func_cos(bmp_real x) {
        return 1.0 / (2 - cos(x));
    }

    bmp_real func_demo(bmp_real x) {
        size_t n = 2;
        bmp_real a = 1.75;
        return (a*a - 1)*pow(a, n)*cos(n*x) / pow(1 - 2 * a*cos(x) + a*a, 1);
    }

    bmp_real func_exp_sin(bmp_real x) {
        return exp(sin(x));
    }

    bmp_real func_exp_cos(bmp_real x) {
        return exp(cos(x));
    }
}

// Сравнение схемы Горнера и метода трапеций epc
TEST_CASE("GornerVsTrapz") {
    std::cout.precision(std::numeric_limits<bmp_real>::max_digits10);
    bmp_real k = 1.0 / 2;
    bmp_real x = -1;
    bmp_real I_base = fdsf::richardson_method(x, 0, k);
    bmp_real I_prec = fdsf::Gorner(x, k);
    std::cout << I_base << std::endl;
    std::cout << I_prec << std::endl;
}

TEST_CASE("ExpConverge") {
    // Для статьи о сверхстепенной сходимости
    //epc::Richardson(func_demo, 0, fdsf::PI, true);
    epc::Richardson(func_cos, 0, 1);
    epc::Richardson(func_exp_cos, 0, 1);
    epc::Richardson(func_exp_sin, 0, 1);
}

TEST_CASE("BOOST_GAMMA") {
    using namespace boost::multiprecision;
    std::setprecision(std::numeric_limits<bmp_real>::max_digits10);
    // Operations at fixed precision and full numeric_limits support:
    cpp_dec_float_100 b = 2;
    std::cout << std::numeric_limits<cpp_dec_float_100>::digits << std::endl;
    // Note that digits10 is the same as digits, since we're base 10! :
    std::cout << std::numeric_limits<cpp_dec_float_100>::digits10 << std::endl;
    // We can use any C++ std lib function, lets print all the digits as well:
    std::cout << std::setprecision(std::numeric_limits<cpp_dec_float_100>::max_digits10)
        << log(b) << std::endl << log(bmp_real(2)) << std::endl; // print log(2)
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
    bmp_real my_pi = boost::math::constants::pi<bmp_real>();
    std::cout << "my_pi = " << my_pi << std::endl;
    std::cout << "4atg1 = " << 4 * boost::multiprecision::atan(bmp_real(1)) << std::endl;
    std::cout << "I1(0) = " << my_pi*my_pi / 12 << std::endl;
}
