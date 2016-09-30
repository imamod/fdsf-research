#pragma once

#include <boost/multiprecision/cpp_dec_float.hpp>

namespace epc{
    using namespace boost::multiprecision;

    typedef cpp_dec_float_50 bmp_real;
    // typedef double bmp_real;

    void checkTrapz(bmp_real a, bmp_real b);
} // excess_power_convergense