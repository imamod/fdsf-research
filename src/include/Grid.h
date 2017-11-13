#pragma once

#include "BasicTypes.h"

class Grid {
    public:
        Grid(size_t N, size_t addNCount = 11);

        // Установить сетку базовых узлов по линейному закону
        void setLinearGrid();

        /**
         * Задает правую линейно-тригонометрическое сетку в базовых узлах слева, плюс m_additionalNCount
         * дополнительных точек между каждой парой базовых узлов.
         * Актуально для для левосторонней аппроксимации целых и полуцелых индексов.
         */
        void setLinearTrigonometricGrid();

        /**
         * Задает правую линейно-тригонометрическое сетку в базовых узлах справа, плюс m_additionalNCount
         * дополнительных точек между каждой парой базовых узлов.
         * Актуально только для полуцелых индексов.
         */
        void setLinearTrigonometricGridRight();

        /**
         * Замена переменных. Исследование сетки с автоматическим выравниванием экстремумов.
         * ksi = 1/y; delta - экстремумы, полученные на сетке с ksi.
         * eta(i) = ksi(i) + 0.5*(ksi(i+1)-ksi(i-1))*tau*(1+sqrt(delta(i-0.5)/delta(i+0.5)))/(1-sqrt(delta(i-0.5)/delta(i+0.5)))
         * tau = 0.5
         */
        void shiftLinTrigGrid(const BmpVector& delta);

        // Получить массив базовых точек
        BmpVector base() const;

        // Получить массив дополнительных точек (вместе с базовыми)
        BmpVector additional() const;

        // Получить x по y
        BmpVector xByY(const BmpVector& y);

    private:
        // Число базовых точек
        size_t m_N_base;
        // Число дополнительных точек между базовыми
        size_t m_additionalNCount;
        // Массив базовых точек
        BmpVector m_base;
        // Массив дополнительных точек
        BmpVector m_additional;

        /* Заполнить сетку дополнительными точками */
        void setAdditionalDots();

        /**
         * В правосторонней аппроксимации для определения распределения узлов интерполяции
         * используется замена переменных ksi = 1/y. Но расчет значения интеграла в точке ведется по переменной y
         * ( точнее говоря по x = ln(exp(y) - 1) ). Поэтому нужны нижеприведенные вспомогательные функции.
         */
        /* Выполнить замену переменных y = 1/ksi */
        void changeGrid(BmpVector& v);
};
