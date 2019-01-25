#pragma once

#include "BasicTypes.h"
#include "Trapz.h"

/* Класс описывающий метод трапеций для функций ФД */
class TrapzFD : public Trapz {
    public:
        TrapzFD(fdsf::SubIntegralFunc func);
        // Вычисление на первой сетке
        virtual double firstGrid(BmpReal initialGrid, const fdsf::Params& params);
        // Экономичное вычисление при следующих сгущениях
        virtual double economicGrid(BmpReal initialGrid, double previousValue, const fdsf::Params& params);

};
