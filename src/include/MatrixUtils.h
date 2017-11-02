#pragma once

#include "BasicTypes.h"
#include "Fdsf.h"
#include <iostream>

class CMatrix {
    public:
        CMatrix(const BmpMatrix& matrix);

        // Заполнение матриц... Переделать
        void fill_matrix(const size_t N_base, BmpVector z,
                         BmpVector y0, BmpVector &B,
                         BmpMatrix &A);

        // Получить аппроксимирующие коэффициенты
        void find_coefficients(BmpMatrix A_inv, BmpVector B,
                               BmpVector &a, BmpVector &b,
                               size_t N);

        // Получить обратную матрицу
        BmpMatrix inverse();

        // Распечатать матрицу
        void print(const BmpMatrix& matrix);

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
