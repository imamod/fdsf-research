#pragma once

#include "BasicTypes.h"
#include "json.hpp"

namespace quad {
    // Вычислить значение ФД индека k для 0 <= x <= x_min
    nlohmann::json calculate(double k, double x);

    // TODO: Вычислить значение интегральной ФД для 0 <= x <= x_min
    BmpReal calculateJmhalf(BmpReal x);
}
