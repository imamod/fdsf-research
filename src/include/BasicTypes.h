#pragma once

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/detail/default_ops.hpp>
#include <boost/multiprecision/number.hpp>
#include <boost/math/special_functions/gamma.hpp>

#include <map>
#include <string>
#include <vector>

//using BmpReal = boost::multiprecision::cpp_dec_float_50;
using bmp_real = boost::multiprecision::cpp_dec_float_50;
using Real = double;

using BmpVector = std::vector<bmp_real>;

namespace fdsf {
    // Желаемая точность расчета
    //const BmpReal epsilon = 1e-45;
    const bmp_real epsilon = 1e-17;

    // Значение pi
    //const BmpReal PI = boost::math::constants::pi<BmpReal>();
    const bmp_real PI = boost::math::constants::pi<bmp_real>();
}
