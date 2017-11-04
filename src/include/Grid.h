#pragma once

#include "BasicTypes.h"

class Grid {
    public:
        Grid(size_t N, size_t addNCount = 11);

        // ”становить сетку базовых узлов по линейному закону
        void setLinearGrid();

        // ”становить сетку базовых узлов линейно-тригонометрическому закону слева
        void setLinearTrigonometricGrid();

        /**
         * «адает правую линейно-тригонометрическое сетку в базовых узлах, плюс 10(?)
         * дополнительных точек между каждой парой базовых узлов. јктуально только дл€ полуцелых индексов
         */
        void setLinearTrigonometricGridRight();

        // ѕолучить массив базовых точек
        BmpVector base() const;

        // ѕолучить массив дополнительных точек (вместе с базовыми)
        BmpVector additional() const;

        // ѕолучить x по y
        BmpVector xByY(const BmpVector& y);

    private:
        // „исло базовых точек
        size_t m_N_base;
        // „исло дополнительных точек между базовыми
        size_t m_additionalNCount;
        // ћассив базовых точек
        BmpVector m_base;
        // ћассив дополнительных точек
        BmpVector m_additional;
};
