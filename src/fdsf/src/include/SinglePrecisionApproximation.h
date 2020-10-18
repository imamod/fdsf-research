#pragma once

#include "BasicTypes.h"

namespace single_precision_formula {
    // Левая аппроксимация
    BmpReal left(BmpReal k, BmpReal x);
    // Правая аппроксимация
    BmpReal right(BmpReal k, BmpReal x);
}
