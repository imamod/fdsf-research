#pragma once

#include "BasicTypes.h"

/* Класс, описывающий сгущение по Ричардсону */
class Richardson {
    public:
        explicit Richardson(double x, std::function<double(double, size_t)> f, size_t initialGrid);

        /* Вычислить значение на сгущающихся сетках */
        void calculate();

        /* Получить вычисленное значение интеграла*/
        const double get() const;

    private:
        // Текущая сетка вычисления
        size_t m_N;
        // Точка в которой вычисляется значение интеграла
        double m_x;
        // Финальное значение интеграла в точке x
        double m_I;
        // Подынтегральная функция
        std::function<double(double, size_t)> m_f;
};
