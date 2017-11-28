#pragma once

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/detail/default_ops.hpp>
#include <boost/multiprecision/number.hpp>
#include <boost/math/special_functions/gamma.hpp>

#include <map>
#include <string>
#include <vector>

#ifdef HIGH_PRECISION
using BmpReal = boost::multiprecision::cpp_dec_float_50;
namespace fdsf {
    // Желаемая точность расчета
    const BmpReal epsilon50 = 1e-45;
}
#else
using BmpReal = double;
#endif

using BmpVector = std::vector<BmpReal>;
using BmpMatrix = std::vector<BmpVector>;

namespace fdsf {
    // Желаемая точность расчета
    const BmpReal epsilon = 1e-17;
    // Значение pi
    const BmpReal PI = boost::math::constants::pi<BmpReal>();
}
