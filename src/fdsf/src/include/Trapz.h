#pragma once

#include "BasicTypes.h"

using FermiFunction = std::function<double(double, double, double)>;

struct FermiDirakFunction {
    // Индекс функции ФД
    double index;
    // Точка x в которой вычисляется значение x
    double x;
    // Представление подынтегральной функции
    std::function<double(double, double, double)> func;
};

/* Класс описывающий метод трапеций для функций ФД */
class TrapzFD {
    public:
        explicit TrapzFD(const FermiDirakFunction& fd, double grid);

        // Вычислить методом трапеций значение интеграла на сетке
        double trapz(double previousGridValue);

    private:
        // Число узлов сетки
        double m_N;
        // Функция ФД
        FermiDirakFunction m_fd;
        // Вычисление на первой сетке
        double firstGrid();
        // Экономичное вычисление при следующих сгущениях
        double economicGrid(double previousValue);
};
