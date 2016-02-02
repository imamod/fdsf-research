// Implementation of matrix_helper.h
//
#include "stdafx.h"
#include "matrix_helper.h"
#include <limits>
#include <iomanip>

double CMatrix::gaus_det(matrix_type::_matrix mass, size_t cnt_str)
{
    double det = 1;
    //прямой ход
    for (size_t i = 0; i < cnt_str; i++)
    {
        for (size_t j = i + 1; j < cnt_str; j++)
        {
            if (mass.at(i).at(i) == 0)
                return 0;
            double b = mass.at(j).at(i) / mass.at(i).at(i);
            for (size_t k = i; k < cnt_str; k++)
                mass.at(j).at(k) = mass.at(j).at(k) - mass.at(i).at(k) * b;
        }
        det *= mass.at(i).at(i);//вычисление определителя
    }
    return det;
}

std::vector < std::vector <double> > inverse(std::vector < std::vector <double> > a) {
    int n = a.size();
    std::vector < std::vector <double> > ans(n, std::vector <double>(n, 0));
    for (int i = 0; i < n; i++){
        ans.at(i).at(i) = 1.0;
    }
    for (int i = 0; i < n; i++){
        int row = i;
        double mx = a.at(i).at(i);
        for (int k = i + 1; k < n; k++){
            if ((abs(a.at(k).at(i)) - mx) > epsilon){
                row = k;
                mx = abs(a.at(k).at(i));
            }
        }
        //if (row != i) {
        //    swap(a[row], a[i]);
        //    swap(ans[row], ans[i]);
        //}
        for (int j = i + 1; j < n; j++){
            double e = a.at(j).at(i) / a.at(i).at(i);
            for (int k = 0; k < n; k++){
                a.at(j).at(k) -= e*a.at(i).at(k);
                ans.at(j).at(k) -= e*ans.at(i).at(k);
            }
        }
    }
    for (int i = n - 1; i >= 0; i--) {
        for (int j = i - 1; j >= 0; j--){
            double e = a.at(j).at(i) / a.at(i).at(i);
            for (int k = 0; k < n; k++){
                a.at(j).at(k) -= e*a.at(i).at(k);
                ans.at(j).at(k) -= e*ans.at(i).at(k);
            }
        }
        for (int j = 0; j < n; j++) {
            ans.at(i).at(j) /= a.at(i).at(i);
            std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 2) << ans.at(i).at(i) << " ";
        }
        std::cout << std::endl;
    }
    return ans;
}

void CMatrix::gaus_inv(matrix_type::_matrix A, size_t size, matrix_type::_matrix &A_inv)
{
    for (size_t i = 0; i < size; i++) {
        matrix_type::_vector ivec;
        for (size_t j = 0; j < size; j++) {
            ivec.push_back(0);
        }
        A_inv.push_back(ivec);
        ivec.clear();
        A_inv.at(i).at(i) = 1.0;
    }

    //прямой ход
    double a, b;
    for (size_t i = 0; i < size; i++)
    {
        a = A.at(i).at(i);
        for (size_t j = i + 1; j < size; j++)
        {
            b = A.at(j).at(i);
            for (size_t k = 0; k < size; k++) {
                A.at(j).at(k) = A.at(j).at(k)*b - A.at(i).at(k) * a;
                A_inv.at(j).at(k) = A_inv.at(j).at(k)*b - A_inv.at(i).at(k) * a;
            }
            std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 2) << A_inv.at(i).at(j) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "\n\n";
    //обратный ход вычисления элементов обратной матрицы
    for (size_t i = 0; i < size; i++)
    {
        for (size_t j = size - 1; j > 0; j--)
        {
            double sum = 0;
            for (size_t k = size - 1; k > j; k--) {
                sum += A.at(j).at(k)*A_inv.at(k).at(i);
            }
            A_inv.at(j).at(i) = (A_inv.at(j).at(i) - sum) / A.at(j).at(j);
            std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 2) << A_inv.at(i).at(j) << " ";
        }
        std::cout << std::endl;
    }

}

void CMatrix::fill_matrix(const int N_base, matrix_type::_vector z, matrix_type::_vector y0,
                          matrix_type::_vector &B, matrix_type::_matrix &A)
{
    for (size_t i = 0; i < 2 * N_base + 1; i++) {
        B.push_back(z.at(i) - 1);
        matrix_type::_vector ivec;
        for (size_t j = 0; j < 2 * N_base + 1; j++) {
            ivec.push_back(0);
        }
        A.push_back(ivec);
        ivec.clear();
    }

    for (size_t i = 0; i < 2 * N_base + 1; i++) {
        for (size_t j = 0; j < 2 * N_base + 1; j++) {
            if (j >= 0 && j <= N_base)
                A.at(i).at(j) = pow(y0.at(i), j + 1 );
            else if (j >= N_base + 1 && j <= 2 * N_base)
                A.at(i).at(j) = -z.at(i)*pow(y0.at(i), j - N_base);
        }
    }
}

void CMatrix::find_coefficients(matrix_type::_matrix A_inv, matrix_type::_vector B,
                                matrix_type::_vector &a, matrix_type::_vector &b, int N)
{
    matrix_type::_vector ksi;
    for (size_t i = 0; i < B.size(); i++) {
        ksi.push_back(0);
        for (size_t j = 0; j < B.size(); j++) {
            ksi.at(i) += A_inv.at(i).at(j)*B.at(j);
            //std::cout << "i = " << i << " j = " << j << ": " << std::setprecision(std::numeric_limits<double>::digits10 + 2) << ksi.at(i) << " " << std::endl;
        }
        //std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 2) << ksi.at(i) << " " << std::endl;
    }

    for (size_t j = 0; j < ksi.size(); j++) {
        if (j >= 0 && j < N + 1) {
            a.push_back(0);
            a.at(j) = ksi.at(j);
        }
        else if (j >= N + 1 && j <= 2 * N) {
            b.push_back(0);
            b.at(j - N - 1) = ksi.at(j);
        }
    }
}

std::ostream& operator << (std::ostream& output, CMatrix& a)
{
    output << a << " ";
    return output;
}
