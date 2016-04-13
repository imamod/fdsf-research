#pragma once

#include "../FDSFInteger.h"
#include <boost\multiprecision\cpp_dec_float.hpp>
#include <iostream>
#include <vector>

//using namespace boost::multiprecision;
using namespace fdsf;

namespace matrix_type{
    
    typedef std::vector<std::vector<cpp_dec_float_50>> _matrix;
    typedef std::vector<cpp_dec_float_50> _vector;
}

std::vector < std::vector <cpp_dec_float_50> > inverse(std::vector < std::vector <cpp_dec_float_50> > a);

class CMatrix {

public:
   // CMatrix();
   // ~CMatrix();

    cpp_dec_float_50 gaus_det(matrix_type::_matrix mass, size_t cnt_str);
    void CMatrix::gaus_inv(matrix_type::_matrix A, size_t size, matrix_type::_matrix &A_inv);
    void fill_matrix(const int N_base, matrix_type::_vector z, matrix_type::_vector y0, matrix_type::_vector &B, matrix_type::_matrix &A);
    void find_coefficients(matrix_type::_matrix A_inv, matrix_type::_vector B, matrix_type::_vector &a, matrix_type::_vector &b, int N);

private:
    friend std::ostream& operator << (std::ostream&, CMatrix& a);
};

void GetApproxomateValues(matrix_type::_vector &a,
    matrix_type::_vector &b,
    matrix_type::_vector &y0,
    matrix_type::_vector &Y,
    matrix_type::_vector &I,
    matrix_type::_vector &z,
    matrix_type::_vector &delta_base,
    matrix_type::_vector &delta_additional, const int N_base);