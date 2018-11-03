#pragma once

#include "BasicTypes.h"
#include "json.hpp"

namespace quad {
    nlohmann::json calculate(double k, double x);
    // TODO: remove
    double euler_maclaurin(double x, double k, int N);
}