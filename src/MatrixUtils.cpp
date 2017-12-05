#include "MatrixUtils.h"
#include "BasicService.h"
#include <math.h>
#include <limits>
#include <iomanip>

MatrixUtils::MatrixUtils(const BmpMatrix& matrix)
    : m_matrix(matrix) {}

/**
 * Сформировать единичную матрицу
 */
BmpMatrix MatrixUtils::eye(const size_t N) {
    BmpMatrix eyeMatrix(N, BmpVector(N, 0));
    for (size_t i = 0; i < N; ++i) {
        eyeMatrix[i][i] = 1.0;
    }
    return eyeMatrix;
}
/**
 * Функция получения аппроксимированный значений z. Использовалось для целых
 * TODO: причесать.
*/
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

/**
 * Получить обратную квадратную матрицу методом Гаусса
 */
BmpMatrix MatrixUtils::inverse() {
    size_t n = m_matrix.size();
    // Формируем единичную матрицу
    BmpMatrix inversedMatrix = eye(n);
    // Прямой ход
    for (size_t i = 0; i < n; ++i) {
        BmpReal mx = m_matrix[i][i];
        // Что делает этот проход???...
        for (size_t k = i + 1; k < n; ++k) {
            // Что за условие???...
            if ((abs(m_matrix[k][i]) - mx) > fdsf::epsilon) {
                mx = abs(m_matrix[k][i]);
            }
        }
        // Приведение матрицы к треугольному виду
        for (size_t j = i + 1; j < n; ++j) {
            // Делим каждый элемент строки на диагональный элемент
            BmpReal e = m_matrix[j][i] / m_matrix[i][i];
            // Пробегаем по столбцу и зануляем все под диагональным элементом
            for (size_t k = 0; k < n; ++k) {
                m_matrix[j][k] -= e*m_matrix[i][k];
                inversedMatrix[j][k] -= e*inversedMatrix[i][k];
            }
        }
    }

    // Обратный ход
    for (int i = n - 1; i >= 0; i--) {
        for (int j = i - 1; j >= 0; j--) {
            BmpReal e = m_matrix[j][i] / m_matrix[i][i];
            for (size_t k = 0; k < n; ++k) {
                m_matrix[j][k] -= e*m_matrix[i][k];
                inversedMatrix[j][k] -= e*inversedMatrix[i][k];
            }
        }
        for (size_t j = 0; j < n; ++j) {
            inversedMatrix[i][j] /= m_matrix[i][i];
        }
    }

    return inversedMatrix;
}
// TODO: Унифицировать. Использовалось для целых
void MatrixUtils::fill_matrix(const size_t N_base, BmpVector z,
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
// TODO: Унифицировать. Использовалось для целых
void MatrixUtils::find_coefficients(BmpMatrix A_inv, BmpVector B,
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
std::ostream& operator << (std::ostream& output, MatrixUtils& a) {
    output << a << " ";
    return output;
}

/*****************************************************************************
 * Solver for right system of half integer *
 *****************************************************************************/

BmpVector solveRightApproximationSystem(BmpReal k, size_t N, const BmpVector& y0, const BmpVector& I_base) {
    const auto baseSize = 2 * N;
    const BmpReal C1 = (k + 1)*k*fdsf::PI*fdsf::PI / 6;
    // Задаем Матрицы А, B и z
    BmpVector B(baseSize, 0);
    BmpVector z(baseSize, 0);
    BmpMatrix A = BmpMatrix(baseSize, BmpVector(baseSize, 0));
    // Вычисляем матрицы А, B и z
    for (int i = 0; i < baseSize; ++i) {
        auto y0_i = y0[i];
        auto underPow = I_base[i] * (k + 1) / pow(y0_i, k + 1);
        z[i] = (pow(underPow, 2 / k) - 1)*y0_i*y0_i*k / (2 * C1);
        B[i] = z[i] - 1;
        for (int j = 0; j < baseSize; ++j) {
#ifdef HIGH_PRECISION
            if (j < N ) {
                A[i][j] = pow(y0_i, -2 * (j + 1));
            } else {
                A[i][j] = -z[i] * pow(y0_i, 2 * ((int)N - j + 1));
            }
#else
            A[i][j] = j < N ? pow(y0_i, -2 * (j + 1))
                : -z[i] * pow(y0_i, 2 * ((int)N - j + 1));
#endif
        }
    }

    // Берем обратную матрицу
    BmpMatrix A_inv = MatrixUtils(A).inverse();

    // Решаем СЛАУ A*E = B
    BmpVector E(baseSize, 0);
    for (size_t i = 0; i < baseSize; ++i) {
        for (size_t j = 0; j < baseSize; ++j) {
            E[i] += A_inv[i][j]*B[j];
        }
    }
    return E;
}

/**
* Получить вектор приближенных значений в точках у для полуцелых индексов k, правая аппроксимация
*/
BmpVector approximateValueRight(const BmpVector& a, const BmpVector& b, const BmpVector& y, BmpReal k) {
    auto N = y.size() / 2;
    const BmpReal C1 = (k + 1)*k*pow(fdsf::PI, 2) / 6;
    BmpVector I;
    for (auto it : y) {
        BmpReal S1 = 1, S2 = 1;
        for (auto n = 1; n < N + 1; ++n) {
            auto y_pow = pow(it, (-2 * n));
            S1 += a[n-1]*y_pow;
            S2 += b[n-1]*y_pow;
        }
        auto F_z_base = S1 / S2;
        I.push_back(pow(((F_z_base * 2 * C1 / (k*it*it)) + 1), (k / 2))*pow(it, k + 1) / (k + 1));
    }
    return I;
}