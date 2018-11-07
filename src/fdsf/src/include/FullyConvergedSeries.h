#pragma once

#include "BasicTypes.h"

namespace fcs {
    // Вычислить значение ФД индекса k для x <= 0
    BmpReal calculate(BmpReal k, BmpReal x);

    // Вычислить значение интегральной ФД для x <= 0
    BmpReal calculateJmhalf(BmpReal x);
}
