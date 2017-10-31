#pragma once

#include "BasicTypes.h"

namespace filesys {
    // Записать данные в файл
    void writeFile(const std::string& filename, const BmpVector& data);
}

// Сверхстепенная сходимость
namespace epc {
    typedef BmpReal(*function)(BmpReal x);

    BmpReal Richardson(function f, BmpReal a, BmpReal b, bool countEvery = false);
} // exponential_convergense
