#include "MatrixUtils.h"
#include <math.h>
#include <limits>
#include <iomanip>

CMatrix::CMatrix(const BmpMatrix& matrix)
    : m_matrix(matrix) {}


void GetApproxomateValues(BmpVector &a,
                          BmpVector &b,
                          BmpVector &y0,
                          BmpVector &Y,
                          BmpVector &I,
                          BmpVector &z,
                          BmpVector &delta_base,
                          BmpVector &delta_additional,
                          const size_t N_base) {
    const size_t baseSize = y0.size();
    BmpVector F_base(baseSize, 0);

    for (size_t j = 0; j < baseSize; j++) {
        BmpReal S1 = 0, S2 = 0;

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
    BmpVector F(addSize, 0);
        

    for (size_t j = 0; j < addSize; j++) {
        BmpReal S1 = 0, S2 = 0;

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

// Получить обратную матрицу
BmpMatrix CMatrix::inverse() {
    size_t n = m_matrix.size();
    BmpMatrix ans(n, BmpVector(n, 0));
    for (size_t i = 0; i < n; i++) {
        ans.at(i).at(i) = 1.0;
    }
    for (size_t i = 0; i < n; i++) {
        size_t row = i;
        BmpReal mx = m_matrix.at(i).at(i);
        for (size_t k = i + 1; k < n; k++) {
            if ((abs(m_matrix.at(k).at(i)) - mx) > fdsf::epsilon) {
                row = k;
                mx = abs(m_matrix.at(k).at(i));
            }
        }

        for (size_t j = i + 1; j < n; j++) {
            BmpReal e = m_matrix.at(j).at(i) / m_matrix.at(i).at(i);
            for (size_t k = 0; k < n; k++){
                m_matrix.at(j).at(k) -= e*m_matrix.at(i).at(k);
                ans.at(j).at(k) -= e*ans.at(i).at(k);
            }
        }
    }
    for (size_t i = n - 1; i >= 0; i--) {
        for (size_t j = i - 1; j >= 0; j--) {
            BmpReal e = m_matrix.at(j).at(i) / m_matrix.at(i).at(i);
            for (size_t k = 0; k < n; k++) {
                m_matrix.at(j).at(k) -= e*m_matrix.at(i).at(k);
                ans.at(j).at(k) -= e*ans.at(i).at(k);
            }
        }
        for (size_t j = 0; j < n; j++) {
            ans.at(i).at(j) /= m_matrix.at(i).at(i);
            //std::cout << std::setprecision(std::numeric_limits<BmpReal>::digits10 + 2) << ans.at(i).at(i) << " ";
        }
        //std::cout << std::endl;
    }
    return ans;
}

void CMatrix::fill_matrix(const size_t N_base, BmpVector z,
                          BmpVector y0,
                          BmpVector &B, BmpMatrix &A) {
    for (size_t i = 0; i < 2 * N_base + 1; i++) {
        B.push_back(z.at(i) - 1);
        BmpVector ivec;
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

void CMatrix::find_coefficients(BmpMatrix A_inv, BmpVector B,
                                BmpVector &a, BmpVector &b, size_t N) {
    BmpVector ksi;
    for (size_t i = 0; i < B.size(); i++) {
        ksi.push_back(0);
        for (size_t j = 0; j < B.size(); j++) {
            ksi.at(i) += A_inv.at(i).at(j)*B.at(j);
            //std::cout << "i = " << i << " j = " << j << ": " << std::setprecision(std::numeric_limits<BmpReal>::digits10 + 2) << ksi.at(i) << " " << std::endl;
        }
        //std::cout << std::setprecision(std::numeric_limits<BmpReal>::digits10 + 2) << ksi.at(i) << " " << std::endl;
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

// Вывод объекта СMatrix
std::ostream& operator << (std::ostream& output, CMatrix& a) {
    output << a << " ";
    return output;
}

// Распечатать матрицу
void CMatrix::print(const BmpMatrix& matrix) {
    for (const auto& row : matrix) {
        for (const auto& item : row) {
            std::cout << std::setw(8) << std::setfill(' ') << item << " ";
        }
        std::cout << std::endl;
    }
}
