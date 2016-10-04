#pragma once

#include <boost\multiprecision\cpp_dec_float.hpp>

namespace newton {
    using namespace boost::multiprecision;

    typedef cpp_dec_float_50 bmp_real;
    //typedef double bmp_real;

    // Метод Ньютона
    bmp_real NewtonsMethod(bmp_real x, bmp_real k);

} //newton
