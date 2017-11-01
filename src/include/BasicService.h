#pragma once

#include "BasicTypes.h"
#include <functional>
#include <iostream>

namespace filesys {
    // Записать данные в файл
    void writeFile(const std::string& filename, const BmpVector& data);
}

// Сверхстепенная сходимость
namespace epc {
    BmpReal Richardson(std::function<BmpReal(const BmpReal&)> f, BmpReal a, BmpReal b, bool countEvery = false);
} // exponential_convergense

namespace fdsf {
    // Метод Ньютона
    BmpReal NewtonsMethod(BmpReal x, BmpReal k);
}
