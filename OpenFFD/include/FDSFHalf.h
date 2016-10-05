#pragma once

#include <vector>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/detail/default_ops.hpp>
#include <boost/multiprecision/number.hpp>
#include <boost/math/special_functions/gamma.hpp>

namespace fdsf {
    using namespace boost::multiprecision;

    typedef cpp_dec_float_50 bmp_real;
    //typedef double bmp_real;

    bmp_real fermi_dirak_half_integer(bmp_real ksi, bmp_real x, bmp_real k, bmp_real a);

    // Formula Euler-Maclaurin 
    bmp_real euler_maclaurin_method(bmp_real x, const bmp_real k, int N, bmp_real& a);

} // fdsf