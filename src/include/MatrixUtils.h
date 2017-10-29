#pragma once

#include "Fdsf.h"
#include <boost\multiprecision\cpp_dec_float.hpp>
#include <iostream>
#include <vector>

//using namespace boost::multiprecision;
using namespace fdsf;

namespace matrix_type{
    
    typedef std::vector<BmpVector> _matrix;
    typedef BmpVector _vector;
}

std::vector < std::vector <bmp_real> > inverse(std::vector < std::vector <bmp_real> > a);

class CMatrix {

public:

    void fill_matrix(const size_t N_base, matrix_type::_vector z,
                     matrix_type::_vector y0, matrix_type::_vector &B,
                     matrix_type::_matrix &A);

    void find_coefficients(matrix_type::_matrix A_inv, matrix_type::_vector B,
                           matrix_type::_vector &a, matrix_type::_vector &b,
                           size_t N);

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
    matrix_type::_vector &delta_additional, const size_t N_base);