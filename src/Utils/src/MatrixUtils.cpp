#include "MatrixUtils.h"
#include "BasicService.h"
#include "Constants.h"
#include <math.h>
#include <limits>
#include <iomanip>

/**
 * Сформировать единичную матрицу
 */
BmpMatrix eye(const size_t N) {
    BmpMatrix eyeMatrix(N, BmpVector(N, 0));
    for (size_t i = 0; i < N; ++i) {
        eyeMatrix[i][i] = 1.0;
    }
    return eyeMatrix;
}

/**
 * Функция получения аппроксимированныx значений I. Использовалось для целых
*/
BmpVector approximateFunctionValueIntegerIndex(const BmpVector& a, const BmpVector& b, const BmpVector& y, size_t N_base) {
    const size_t baseSize = y.size();
    BmpVector I(baseSize, 0);

    for (size_t j = 0; j < baseSize; ++j) {
        BmpReal S1 = 1, S2 = 1;

        for (size_t n = 0; n < N_base + 1; ++n) {
            S1 = S1 + a[n]*pow(y[j], n + 1);
        }
        for (size_t m = 0; m < N_base; ++m) {
            S2 = S2 + b[m]*pow(y[j], m + 1);
        }
 
         I[j] = S1 / S2;
    }
    return I;
}

/**
* Получить погрешноть
*/
BmpVector relativeApproximationError(const BmpVector& precisionF, const BmpVector& approximateF) {
    if (precisionF.size() != approximateF.size()) {
        throw std::domain_error("Vectors must be the same length.");
    }
    const size_t n = precisionF.size();
    BmpVector delta(n, 0);
    for (auto i = 0; i < n; ++i) {
        delta[i] = approximateF[i] / precisionF[i] - 1;
    }
    return delta;
}

/**
 * Получить обратную квадратную матрицу методом Гаусса
 */
BmpMatrix inverse(const BmpMatrix& A) {
    BmpMatrix matrix(A);
    size_t n = matrix.size();
    // Формируем единичную матрицу
    BmpMatrix inversedMatrix = eye(n);
    // Прямой ход
    for (size_t i = 0; i < n; ++i) {
        BmpReal mx = matrix[i][i];
        // Что делает этот проход???...
        for (size_t k = i + 1; k < n; ++k) {
            // Что за условие???...
            if ((abs(matrix[k][i]) - mx) > fdsf::epsilon) {
                mx = abs(matrix[k][i]);
            }
        }
        // Приведение матрицы к треугольному виду
        for (size_t j = i + 1; j < n; ++j) {
            // Делим каждый элемент строки на диагональный элемент
            BmpReal e = matrix[j][i] / matrix[i][i];
            // Пробегаем по столбцу и зануляем все под диагональным элементом
            for (size_t k = 0; k < n; ++k) {
                matrix[j][k] -= e*matrix[i][k];
                inversedMatrix[j][k] -= e*inversedMatrix[i][k];
            }
        }
    }

    // Обратный ход
    for (int i = n - 1; i >= 0; i--) {
        for (int j = i - 1; j >= 0; j--) {
            BmpReal e = matrix[j][i] / matrix[i][i];
            for (size_t k = 0; k < n; ++k) {
                matrix[j][k] -= e*matrix[i][k];
                inversedMatrix[j][k] -= e*inversedMatrix[i][k];
            }
        }
        for (size_t j = 0; j < n; ++j) {
            inversedMatrix[i][j] /= matrix[i][i];
        }
    }

    return inversedMatrix;
}

/**
 * Решить систему для аппроксимации функций целых индексов.
 */
BmpVector solveSystemForIntegerIndex(const BmpVector& z, const BmpVector& y0, size_t N) {
    size_t baseSize = 2 * N + 1;
    // Задаем Матрицы А, B и z
    BmpVector B(baseSize, 0);
    BmpMatrix A = BmpMatrix(baseSize, BmpVector(baseSize, 0));

    for (size_t i = 0; i < baseSize; ++i) {
        B[i] = z[i] - 1;
        for (size_t j = 0; j < baseSize; ++j) {
            if (j < N + 1) {
                A[i][j] = pow(y0[i], j + 1);
            } else {
                A[i][j] = -z[i]*pow(y0[i], j - N);
            }
        }
    }

    BmpMatrix A_inv = inverse(A);
    BmpVector ksi(baseSize, 0);
    for (size_t i = 0; i < baseSize; ++i) {
        for (size_t j = 0; j < baseSize; ++j) {
            ksi[i] += A_inv[i][j]*B[j];
            //std::cout << "i = " << i << " j = " << j << ": " << std::setprecision(std::numeric_limits<BmpReal>::digits10 + 2) << ksi.at(i) << " " << std::endl;
        }
        //std::cout << std::setprecision(std::numeric_limits<BmpReal>::digits10 + 2) << ksi.at(i) << " " << std::endl;
    }
    return ksi;
    //for (size_t j = 0; j < ksi.size(); j++) {
    //    if (j >= 0 && j < N + 1) {
    //        a.push_back(0);
    //        a.at(j) = ksi.at(j);
    //    }
    //    else if (j >= N + 1 && j <= 2 * N) {
    //        b.push_back(0);
    //        b.at(j - N - 1) = ksi.at(j);
    //    }
    //}
}

/*****************************************************************************
 * Solver for right system of half integer *
 *****************************************************************************/

BmpVector solveRightApproximationSystem(BmpReal k, size_t N, const BmpVector& y0, const BmpVector& I_base) {
    const auto baseSize = 2 * N;
    const BmpReal C1 = (k + 1)*k*pi()*pi() / 6;
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
    BmpMatrix A_inv = inverse(A);

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
    const BmpReal C1 = (k + 1)*k*pow(pi(), 2) / 6;
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