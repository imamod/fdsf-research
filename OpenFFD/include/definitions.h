#pragma once

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/detail/default_ops.hpp>
#include <boost/multiprecision/number.hpp>
#include <boost/math/special_functions/gamma.hpp>

namespace fdsf {
    using namespace boost::multiprecision;

    typedef cpp_dec_float_50 bmp_real;
    //typedef double bmp_real;

    // Желаемая точность расчета
    const bmp_real epsilon = 1e-45;
    //const bmp_real epsilon = 1e-17;

    // Значение pi
    const bmp_real PI = boost::math::constants::pi<bmp_real>();

} // fdsf