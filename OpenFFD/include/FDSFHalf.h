#pragma once

#include "definitions.h"
#include <vector>

namespace fdsf {

    // for k = -3/2
    bmp_real fermi_dirak_m3half(bmp_real ksi, bmp_real x, bmp_real k, bmp_real a);
    
    // for others half-integer k
    bmp_real fermi_dirak_half_integer(bmp_real ksi, bmp_real x, bmp_real k, bmp_real a);

    // Formula Euler-Maclaurin 
    bmp_real euler_maclaurin_method(bmp_real x, const bmp_real k, int N, bmp_real& a);

} // namespace fdsf