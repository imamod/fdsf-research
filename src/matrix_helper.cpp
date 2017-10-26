// Implementation of matrix_helper.h
//
#include "matrix_helper.h"
#include <math.h>
#include <limits>
#include <iomanip>

void GetApproxomateValues(matrix_type::_vector &a,
                          matrix_type::_vector &b,
                          matrix_type::_vector &y0,
                          matrix_type::_vector &Y,
                          matrix_type::_vector &I,
                          matrix_type::_vector &z,
                          matrix_type::_vector &delta_base,
                          matrix_type::_vector &delta_additional, 
                          const size_t N_base)
{
    const size_t baseSize = y0.size();
    matrix_type::_vector F_base(baseSize, 0);

    for (size_t j = 0; j < baseSize; j++) {
        bmp_real S1 = 0, S2 = 0;

        for (size_t n = 0; n < N_base + 1; n++) {
            S1 = S1 + a.at(n)*pow(y0.at(j), n + 1);
        }
        for (size_t m = 0; m < N_base; m++) {
                S2 = S2 + b.at(m)*pow(y0.at(j), m + 1);
        }
 
         F_base.at(j) = (1 + S1) / (1 + S2);
         delta_base.at(j) = (F_base.at(j) / z.at(j) - 1);
    }
    //-------------------------------------- - 
    //Äîáàâèì âñïîìîãàòåëüíóþ ñåòêó 
    const size_t addSize = Y.size();
    matrix_type::_vector F(addSize, 0);
        

    for (size_t j = 0; j < addSize; j++) {
        bmp_real S1 = 0, S2 = 0;

        for (size_t n = 0; n < N_base + 1; n++) {
            S1 = S1 + a.at(n)*pow(Y.at(j), n + 1);
            
        }
        for (size_t m = 0; m < N_base; m++) {
            S2 = S2 + b.at(m)*pow(Y.at(j), m + 1);
        }

        F.at(j) = (1 + S1) / (1 + S2);
        delta_additional.at(j) = (F.at(j) / I.at(j) - 1);
    }
}

std::vector < std::vector <bmp_real> > inverse(std::vector < std::vector <bmp_real> > a) 
{
    size_t n = a.size();
    std::vector < std::vector <bmp_real> > ans(n, std::vector <bmp_real>(n, 0));
    for (size_t i = 0; i < n; i++){
        ans.at(i).at(i) = 1.0;
    }
    for (auto i = 0; i < n; i++){
        int row = i;
        bmp_real mx = a.at(i).at(i);
        for (auto k = i + 1; k < n; k++){
            if ((abs(a.at(k).at(i)) - mx) > epsilon){
                row = k;
                mx = abs(a.at(k).at(i));
            }
        }

        for (auto j = i + 1; j < n; j++){
            bmp_real e = a.at(j).at(i) / a.at(i).at(i);
            for (auto k = 0; k < n; k++){
                a.at(j).at(k) -= e*a.at(i).at(k);
                ans.at(j).at(k) -= e*ans.at(i).at(k);
            }
        }
    }
    for (auto i = n - 1; i >= 0; i--) {
        for (auto j = i - 1; j >= 0; j--){
            bmp_real e = a.at(j).at(i) / a.at(i).at(i);
            for (auto k = 0; k < n; k++){
                a.at(j).at(k) -= e*a.at(i).at(k);
                ans.at(j).at(k) -= e*ans.at(i).at(k);
            }
        }
        for (auto j = 0; j < n; j++) {
            ans.at(i).at(j) /= a.at(i).at(i);
            //std::cout << std::setprecision(std::numeric_limits<bmp_real>::digits10 + 2) << ans.at(i).at(i) << " ";
        }
        //std::cout << std::endl;
    }
    return ans;
}

void CMatrix::fill_matrix(const size_t N_base, matrix_type::_vector z,
                          matrix_type::_vector y0,
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
                                matrix_type::_vector &a, matrix_type::_vector &b, size_t N)
{
    matrix_type::_vector ksi;
    for (size_t i = 0; i < B.size(); i++) {
        ksi.push_back(0);
        for (size_t j = 0; j < B.size(); j++) {
            ksi.at(i) += A_inv.at(i).at(j)*B.at(j);
            //std::cout << "i = " << i << " j = " << j << ": " << std::setprecision(std::numeric_limits<bmp_real>::digits10 + 2) << ksi.at(i) << " " << std::endl;
        }
        //std::cout << std::setprecision(std::numeric_limits<bmp_real>::digits10 + 2) << ksi.at(i) << " " << std::endl;
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
