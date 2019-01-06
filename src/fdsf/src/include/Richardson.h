#pragma once

#include "BasicTypes.h"
#include "Trapz.h"

#include <memory>

/* Класс, описывающий сгущение по Ричардсону */
class Richardson {
    public:
        explicit Richardson(std::shared_ptr<FermiDirakFunction>& f, size_t initialGrid);

        /* Вычислить значение на сгущающихся сетках */
        void calculate();

        /* Получить вычисленное значение интеграла*/
        const double get() const;

    private:
        // Текущая сетка вычисления
        size_t m_N;
        // Финальное значение интеграла в точке x
        double m_I;
        // Подынтегральная функция ФД
        std::shared_ptr<FermiDirakFunction> m_fd;
};
