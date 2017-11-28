#pragma once

#include "BasicTypes.h"
#include "Fdsf.h"
#include <iostream>

class CMatrix {
    public:
        CMatrix() {}
        CMatrix(const BmpMatrix& matrix);

        // Заполнение матриц... Переделать
        void fill_matrix(const size_t N_base, BmpVector z,
                         BmpVector y0, BmpVector &B,
                         BmpMatrix &A);

        // Получить аппроксимирующие коэффициенты
        void find_coefficients(BmpMatrix A_inv, BmpVector B,
                               BmpVector &a, BmpVector &b,
                               size_t N);

        /**
         * Сформировать единичную матрицу
         */
        BmpMatrix eye(const size_t N);

        // Получить обратную матрицу методом Гаусса
        BmpMatrix inverse();

    private:
        BmpMatrix m_matrix;
        // Оператор потокового вывода
        friend std::ostream& operator << (std::ostream&, CMatrix& a);
};

void GetApproxomateValues(BmpVector &a,
    BmpVector &b,
    BmpVector &y0,
    BmpVector &Y,
    BmpVector &I,
    BmpVector &z,
    BmpVector &delta_base,
    BmpVector &delta_additional, const size_t N_base);

/**
 * Решить систему для правой аппроксимации
 */
void solveRightApproximationSystem(BmpReal k, size_t N, const BmpVector& y0, const BmpVector& I_base);
