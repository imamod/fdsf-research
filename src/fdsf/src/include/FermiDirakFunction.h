#pragma once

#include "BasicTypes.h"
#include "Trapz.h"

#include <memory>

// Класс, описывающий функцию ФД
class FermiDirakFunction {
    public:
        explicit FermiDirakFunction(double x, double index);

        // Вычислить значение ФД в точке x
        BmpReal calculate(BmpReal grid, BmpReal previousGridValue = 0);

        // Возвращает индекс вычисляемой функции
        BmpReal index() const;
    private:
        // Аргументы ФД
        fdsf::Params m_params;

        BmpReal trapz(const std::shared_ptr<Trapz>& trapz, BmpReal grid, BmpReal previousGridValue);
};
