#pragma once

#include "BasicTypes.h"

// ТODO: class series АБК
/**
* Класс расчета ряда log(n)/(n^2), необходим для вычисления коээффициентов при x>>1
*/
class SeriesLogNDivSqrN {
    public:
        SeriesLogNDivSqrN();
        SeriesLogNDivSqrN(size_t optimal_N);

        // Получение значения
        BmpReal get() const;
        // Число членов, по которому обрезают сумму
        size_t upperBound() const;
        // Рассчет суммы 2<=n<=N-1
        BmpReal limitedSum() const;

    private:
        // Оптимальное N
        size_t m_N{ 300 };
        // Значение ряда
        BmpReal m_seriesSum{ 0 };

        // Вычисление суммы ряда
        BmpReal sum() const;
        // Расчет первой поправки
        BmpReal first() const;
        // Расчет второй поправки
        BmpReal second() const;
        // Интегральная часть
        BmpReal integralPart() const;
};
