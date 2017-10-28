#pragma once

#include "BasicTypes.h"

// TODO: remove this file
namespace epc{
    using namespace fdsf;
    typedef bmp_real(*function)(bmp_real x);

    bmp_real Richardson(function f, bmp_real a, bmp_real b, bool countEvery = false);
} // excess_power_convergense