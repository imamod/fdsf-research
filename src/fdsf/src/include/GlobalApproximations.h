#pragma once

#include "BasicTypes.h"

enum GlobalFive {
    ONE_COEFFICIENT = 1,
    LOW_TEMP,
    BEST_PREC
};

namespace global_formula {
    /**
     * Вычисляет функцию ФД по глобальным 5-членным формулам
     * @param точка x, в которой нужно вычислить значение функции ФД
     * @param индекс функции ФД
     * @param тип глобальной формулы
     * @return вычисленное значение функции ФД
     */
    BmpReal calculateFive(BmpReal x, BmpReal k, uint8_t approxType);
}
