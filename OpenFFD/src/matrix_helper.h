#pragma once

#include "mpreal.h"
#include <vector>

// Желаемая точность расчета
const int digits = 30;

namespace matrix_type{
    typedef std::vector<std::vector<mpfr::mpreal>> _matrix;
    typedef std::vector<mpfr::mpreal> _vector;
}

void Zeidel(matrix_type::_vector &b,
            matrix_type::_matrix &a,
            matrix_type::_vector &x,
            const int n);

std::vector < std::vector <mpfr::mpreal> > inverse(std::vector < std::vector <mpfr::mpreal> > a);

class CMatrix {

public:
   // CMatrix();
   // ~CMatrix();

    mpfr::mpreal gaus_det(matrix_type::_matrix mass, size_t cnt_str);
    void fill_matrix(const mpfr::mpreal N_base, matrix_type::_vector z, matrix_type::_vector y0, matrix_type::_vector &B, matrix_type::_matrix &A);
    void find_coefficients(matrix_type::_matrix A_inv, matrix_type::_vector B, matrix_type::_vector &a, matrix_type::_vector &b, mpfr::mpreal N);

private:
    friend std::ostream& operator << (std::ostream&, CMatrix& a);
};

void GetApproxomateValues(matrix_type::_vector &a,
                          matrix_type::_vector& b,
                          matrix_type::_vector y0,
                          matrix_type::_vector &Y,
                          matrix_type::_vector &I,
                          matrix_type::_vector &z, 
                          matrix_type::_vector &delta_base,
                          matrix_type::_vector &delta_additional, const mpfr::mpreal N_base);