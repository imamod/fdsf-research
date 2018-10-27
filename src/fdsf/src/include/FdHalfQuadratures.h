#pragma once

#include "BasicTypes.h"
#include "json.hpp"

namespace fd {
    const std::string X = "x";
    const std::string I = "I";
    const std::string K = "k";
    const std::string N_MAX = "N";
    const std::string RESULT = "result";
}

namespace quad {
    nlohmann::json calculate(double k, double x);
}