#pragma once

#include "BasicTypes.h"

namespace asympt_series {

    // Cтруктура хранения x_min и N_max (см. препринт2, табл 4)
    struct Limits {
        BmpReal x_min;
        int N_max;
    };

    // Получение x_min и N_max для конкретного индекса
    Limits limits(BmpReal k);

    // Вычислить значение ФД для x >= x_min
    BmpReal calculate(BmpReal k, BmpReal x);
}

        // Значение коэффициентов А асимптотического ряда
        /*const BmpVector m_A = {
            pow(pi(), 2) / 6.0,
            pow(pi(), 4) / 90.0,
            pow(pi(), 6) / 945.0,
            pow(pi(), 8) / 9450.0,
            pow(pi(), 10) / 93555.0,
            691.0 * pow(pi(), 12) / 638512875.0
        };*/

