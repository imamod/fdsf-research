#pragma once

#include "BasicTypes.h"

namespace print {
    // Распечатать вектор
    void vector(const BmpVector& v, bool printInColumn = false);
    // Распечатать матрицу
    void matrix(const BmpMatrix& m);
    // Распечатать мапу (удобно для пар x, f(x))
    void map(const std::map<BmpReal, BmpReal>& m);
}

/* Установить расширенный потоковый вывод на 16 знаков */
void setPreciseOutput();
