#pragma once

#include "BasicTypes.h"
#include "Trapz.h"
#include "FermiDirakFunction.h"

#include <memory>

// Структура, описывающая результат выполнения метода Ричардсона
struct RichardsonResult {
    // Финальная сетка вычисления
    double N;
    // Финальное значение интеграла в точке x
    double I;
};

/* Класс, описывающий сгущение по Ричардсону */
class Richardson {
    public:
        explicit Richardson(BmpReal initialGrid, const FermiDirakFunction& fd);

        /* Вычислить значение на сгущающихся сетках */
        RichardsonResult calculate();

    private:
        // Результат выполнения процедуры Ричардсона
        RichardsonResult m_result;
        // Подынтегральная функция ФД
        std::shared_ptr<FermiDirakFunction> m_func;
};
