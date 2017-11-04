#pragma once

#include "BasicTypes.h"
#include <functional>

namespace test {
    void printVector(const BmpVector& v);
}

namespace filesys {
    // Записать данные в файл
    void writeFile(const std::string& filename, const BmpVector& data);
}

// Сверхстепенная сходимость
namespace epc {
    BmpReal Richardson(std::function<BmpReal(const BmpReal&)> f, BmpReal a, BmpReal b, bool countEvery = false);
}

namespace fdsf {
    // Метод Ньютона
    BmpReal NewtonsMethod(BmpReal x, BmpReal k);
}
