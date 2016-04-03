#pragma once

#include <iostream>
#include <vector>

// Желаемая точность расчета
const double epsilon = 1e-17;

namespace matrix_type{
    typedef std::vector<std::vector<double>> _matrix;
    typedef std::vector<double> _vector;
}

std::vector < std::vector <double> > inverse(std::vector < std::vector <double> > a);

class CMatrix {

public:
   // CMatrix();
   // ~CMatrix();

    double gaus_det(matrix_type::_matrix mass, size_t cnt_str);
    void CMatrix::gaus_inv(matrix_type::_matrix A, size_t size, matrix_type::_matrix &A_inv);
    void fill_matrix(const int N_base, matrix_type::_vector z, matrix_type::_vector y0, matrix_type::_vector &B, matrix_type::_matrix &A);
    void find_coefficients(matrix_type::_matrix A_inv, matrix_type::_vector B, matrix_type::_vector &a, matrix_type::_vector &b, int N);

private:
    friend std::ostream& operator << (std::ostream&, CMatrix& a);
};
