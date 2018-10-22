#pragma once

#include "BasicTypes.h"
#include <functional>

// Сверхстепенная сходимость
namespace epc {
    BmpReal Richardson(std::function<BmpReal(const BmpReal&)> f, BmpReal a, BmpReal b, bool countEvery = false);
}

namespace fdsf {
    // Метод Ньютона
    BmpReal NewtonsMethod(BmpReal x, BmpReal k);
}
