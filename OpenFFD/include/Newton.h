#pragma once

#include <boost\multiprecision\cpp_dec_float.hpp>

namespace fdsf {

    typedef double bmp_real;
    typedef bmp_real(*function)(bmp_real x, bmp_real a, bmp_real k);

    // Метод Ньютона
    bmp_real NewtonsMethod(bmp_real x, bmp_real k);

} // fdsf