// Implementation of matrix_helper.h
//
#include "matrix_helper.h"
#include <math.h>
#include <iostream>

// Условие окончания
bool converge(matrix_type::_vector &xk, matrix_type::_vector &xkp, const int n)
{
    mpfr::mpreal norm = 0;
    for (int i = 0; i < n; i++)
    {
        norm += (xk[i] - xkp[i])*(xk[i] - xkp[i]);
    }
    if (sqrt(norm) >= epsilon)
        return false;
    return true;
}

/*
Ход метода, где:
a[n][n] - Матрица коэффициентов
x[n], p[n] - Текущее и предыдущее решения
b[n] - Столбец правых частей
Все перечисленные массивы вещественные и
должны быть определены в основной программе,
также в массив x[n] следует поместить начальное
приближение столбца решений (например, все нули)
*/
void Zeidel(matrix_type::_vector &b, 
            matrix_type::_matrix &a, 
            matrix_type::_vector &x, 
            const int n)
{
    matrix_type::_vector p = {0,0,0,0,0,0,0,0,0};
    do
    {
        for (int i = 0; i < n; i++)
            p[i] = x[i];

        for (int i = 0; i < n; i++)
        {
            mpfr::mpreal var = 0;
            for (int j = 0; j < i; j++)
                var += (a[i][j] * x[j]);
            for (int j = i + 1; j < n; j++)
                var += (a[i][j] * p[j]);
            x[i] = (b[i] - var) / a[i][i];
            std::cout << x[i] << std::endl;
        }
        std::cout << "-----DEBUG-----" << std::endl;
    } while (!converge(x, p, n));
}

void GetApproxomateValues(matrix_type::_vector &a, 
                          matrix_type::_vector&b, 
                          matrix_type::_vector y0, 
                          matrix_type::_vector &Y, 
                          matrix_type::_vector &I, 
                          matrix_type::_vector &z, 
                          matrix_type::_vector &delta_base, 
                          matrix_type::_vector &delta_additional, const mpfr::mpreal N_base)
{
    const size_t baseSize = y0.size();
    matrix_type::_vector F_base (baseSize,0);
    mpfr::mpreal one = 1.0;
    
    for (int j = 0; j < baseSize; j++) {
        mpfr::mpreal S1 = 0, S2 = 0;
//#if 0
        for (int n = 0; n < N_base+1; n++) {
            S1 = S1 + a.at(n)*mpfr::pow(y0.at(j), n+1);
        }
        for (int m = 0; m < N_base; m++) {
            S2 = S2 + b.at(m)*mpfr::pow(y0.at(j), m+1);
        }
//#endif
#if 0 // Попробуем схитрить и просуммировать с конца
        for (int n = N_base-1; n >=0; n--) {
            S1 = S1 + a.at(n)*mpfr::pow(y0.at(baseSize - j - 1), n + 1);
        }
        for (int m = N_base - 2; m >=0 ; m--) {
            S2 = S2 + b.at(m)*mpfr::pow(y0.at(baseSize - j - 1), m + 1);
        }
#endif
        F_base.at(j) = (1 + S1) / (1 + S2);
        delta_base.at(j) = (F_base.at(j) / z.at(j) - 1);
    }
    //-------------------------------------- -
    //Добавим вспомогательную сетку
    const size_t addSize = Y.size();
    matrix_type::_vector F(addSize, 0);
    //delta_additional.push_back(0);

    for (int j = 0; j < addSize; j++) {
        mpfr::mpreal S1 = 0, S2 = 0;
//#if 0
        for (int n = 0; n < N_base + 1 ; n++) {
            S1 = S1 + a.at(n)*mpfr::pow(Y.at(j), n+1);
        }
        for (int m = 0; m < N_base; m++) {
            S2 = S2 + b.at(m)*mpfr::pow(Y.at(j), m+1);
        }
//#endif
#if 0 //То же самое и здесь
        for (int n = N_base-1; n >= 0; n--) {
            S1 = S1 + a.at(n)*mpfr::pow(Y.at(addSize - j-1), n + 1);
        }
        for (int m = N_base -2; m >= 0; m--) {
            S2 = S2 + b.at(m)*mpfr::pow(Y.at(addSize - j-1), m + 1);
        }
#endif
        F.at(j) = (1 + S1) / (1 + S2);
        delta_additional.at(j) = (F.at(j) / I.at(j) - 1);
    }
}

// вычисление определителя
mpfr::mpreal CMatrix::gaus_det(matrix_type::_matrix mass, size_t cnt_str)
{
    mpfr::mpreal det = 1.0;
    //прямой ход
    for (size_t i = 0; i < cnt_str; i++)
    {
        for (size_t j = i + 1; j < cnt_str; j++)
        {
            if (mass.at(i).at(i) == 0)
                return 0;
            mpfr::mpreal b = mass.at(j).at(i) / mass.at(i).at(i);
            for (size_t k = i; k < cnt_str; k++)
                mass.at(j).at(k) = mass.at(j).at(k) - mass.at(i).at(k) * b;
        }
        det *= mass.at(i).at(i);//вычисление определителя
    }
    return det;
}

// обратная матрица(gauss)
std::vector < std::vector <mpfr::mpreal> > inverse(std::vector < std::vector <mpfr::mpreal> > a) 
{
    int n = a.size();
    std::vector < std::vector <mpfr::mpreal> > ans(n, std::vector <mpfr::mpreal>(n, 0));
    
    std::cout.precision(digits);
    // Initial
    for (int i = 0; i < n; i++){
        ans.at(i).at(i) = 1.0;
    }

    // to triangle view
    for (int i = 0; i < n; i++){
        int row = i;
        mpfr::mpreal mx = a.at(i).at(i);
        // условие выбора главного элемента здесь нам не нужно, потому что и так они на диагонали
       /* for (int k = i + 1; k < n; k++){ 
            if ((abs(a.at(k).at(i)) - mx) > epsilon){
                row = k;
                mx = abs(a.at(k).at(i));
                std::cout << "in condition mx" << std::endl;
            }
        }*/

        for (int j = i + 1; j < n; j++){
            mpfr::mpreal e = a.at(j).at(i) / a.at(i).at(i);
            for (int k = 0; k < n; k++){
                a.at(j).at(k) -= e*a.at(i).at(k);
                ans.at(j).at(k) -= e*ans.at(i).at(k);
            }
        }
    }
    // backword
    for (int i = n - 1; i >= 0; i--) {
        for (int j = i - 1; j >= 0; j--){
            mpfr::mpreal e = a.at(j).at(i) / a.at(i).at(i);
            for (int k = 0; k < n; k++){
                a.at(j).at(k) -= e*a.at(i).at(k);
                ans.at(j).at(k) -= e*ans.at(i).at(k);
            }
        }

        for (int j = 0; j < n; j++) {
            ans.at(i).at(j) /= a.at(i).at(i);
            std::cout << ans.at(i).at(j) << " ";
        }
        std::cout << std::endl;
    }

    return ans;
}

// Заполнение матрицы коэффициентов
void CMatrix::fill_matrix(const mpfr::mpreal N_base, matrix_type::_vector z, matrix_type::_vector y0,
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
            if (j >= 0 && j <= N_base) {
                A.at(i).at(j) = mpfr::pow(y0.at(i), j + 1);
            } else if (j >= N_base + 1 && j <= 2 * N_base){
                A.at(i).at(j) = -z.at(i)*mpfr::pow(y0.at(i), j - N_base);
            }
        }
    }
}

// выделение коэффициентов
void CMatrix::find_coefficients(matrix_type::_matrix A_inv, matrix_type::_vector B,
                                matrix_type::_vector &a, matrix_type::_vector &b, mpfr::mpreal N)
{
    matrix_type::_vector ksi;
    std::cout.precision(digits);
    for (size_t i = 0; i < B.size(); i++) {
        ksi.push_back(0);
        for (size_t j = 0; j < B.size(); j++) {
            ksi.at(i) += A_inv.at(i).at(j)*B.at(j);
            //std::cout << "i = " << i << " j = " << j << ": " << ksi.at(i) << " " << std::endl;
        }
        std::cout << ksi.at(i) << " " << std::endl;
    }

    for (size_t j = 0; j < ksi.size(); j++) {
        if (j >= 0 && j < N + 1) {
            a.push_back(0);
            a.at(j) = ksi.at(j);
        }
        else if (j >= N + 1 && j <= 2 * N) {
            b.push_back(0);
            b.at(j - (int)N - 1) = ksi.at(j);
        }
    }
}

std::ostream& operator << (std::ostream& output, CMatrix& a)
{
    output << a << " ";
    return output;
}
