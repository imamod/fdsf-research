#pragma once

#include "BasicTypes.h"

namespace fdsf {
    using SubIntegralFunc = std::function<double(double, double, double)>;

    struct Params {
        // Индекс функции ФД
        double index;
        // Точка x в которой вычисляется значение x
        double x;
    };
}


class Trapz {
    public:
        explicit Trapz(fdsf::SubIntegralFunc func) 
            : m_func(func) {}
        virtual ~Trapz() {}

        // Вычисление на первой сетке
        virtual double firstGrid(BmpReal initialGrid, const fdsf::Params& params) = 0;
        // Экономичное вычисление при следующих сгущениях
        virtual double economicGrid(BmpReal initialGrid, double previousValue, const fdsf::Params& params) = 0;

        inline fdsf::SubIntegralFunc Trapz::subIntegralFunction() const {
            return m_func;
        }

    private:
        // Колбек подынтегральной функции
        fdsf::SubIntegralFunc m_func;
};
