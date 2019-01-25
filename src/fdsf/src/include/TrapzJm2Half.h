#pragma once

#include "BasicTypes.h"
#include "Trapz.h"

/* Класс описывающий метод трапеций для интегральной функции ФД */
class TrapzJm2Half : public Trapz {
    public:
        TrapzJm2Half(fdsf::SubIntegralFunc func);
        /* Прямая реализация метода трапеций (формула (?) в препринте 2)*/
        virtual double firstGrid(BmpReal initialGrid, const fdsf::Params& params);
        /* Экономичная реализация метода трапеций */
        virtual double economicGrid(BmpReal initialGrid, double previousGridValue, const fdsf::Params& params);

    private:
        /* Веса функции */
        double weights(int n, int m);
};
