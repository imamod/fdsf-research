#pragma once

#include "BasicTypes.h"
#include "json.hpp"

namespace quad {
    nlohmann::json calculate(double k, double x);
}
