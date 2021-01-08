#pragma once

#include "BasicTypes.h"

namespace fdsf {


    // Сгущение по Ричардсону результата работы функции gauss_christoffel_method
    BmpReal richardson_method(BmpReal x, BmpReal k, BmpReal t = 0);



} // namespace fdsf
