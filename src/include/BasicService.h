#pragma once

#include "BasicTypes.h"
#include <functional>

namespace filesys {
    // Считать вектор данных из файла
    BmpVector readFile(const std::string& filename);
    // Записать данные в файл
    void writeFile(const std::string& filename, const BmpVector& data);
    // Создать папку
    std::string createDirectory(BmpReal k, size_t n, const std::string& prefix = "");
}

// Сверхстепенная сходимость
namespace epc {
    BmpReal Richardson(std::function<BmpReal(const BmpReal&)> f, BmpReal a, BmpReal b, bool countEvery = false);
}

namespace fdsf {
    // Метод Ньютона
    BmpReal NewtonsMethod(BmpReal x, BmpReal k);
}
