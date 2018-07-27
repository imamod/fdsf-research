#pragma once

#include "BasicService.h"
#include "Constants.h"

class AsymptoticSeries {
    public:
        AsymptoticSeries(BmpReal k, BmpReal x);

        // Получить вычисленное значение суммы ряда в точке Х
        BmpReal get() const;

    private:
        // Оптимальное N
        size_t m_N{ 6 };
        // Значение коэффициентов А асимптотического ряда
        const BmpVector m_A = {
            pow(pi(), 2) / 6.0,
            pow(pi(), 4) / 90.0,
            pow(pi(), 6) / 945.0,
            pow(pi(), 8) / 9450.0,
            pow(pi(), 10) / 93555.0,
            691.0 * pow(pi(), 12) / 638512875.0
        };
        // Значение ряда
        BmpReal m_seriesSum{ 0 };

        // Рассчитать коэффициенты С
        BmpVector getC(BmpReal k, BmpReal x);
        // Просуммировать
        BmpReal sum(BmpReal k, BmpReal x);

};
