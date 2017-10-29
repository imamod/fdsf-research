#include "TestCommon.h"
#include "Fdsf.h"
#include <unordered_map>

namespace {
    // Для целых
    BmpVector computeIntegral(const BmpVector& x, size_t k) {
        // std::cout << "      I1(0) : " << fdsf::PI*fdsf::PI / 12 << std::endl;
        /**
         * Соответствие значений t для каждого k, подобрано экспериментально
         * bmp_real t = fdsf::get_T_max(X.at(i), k);
         */
        std::unordered_map<size_t, bmp_real> K_T_MAP = {
            { 1, 60 },
            { 2, 75 },
            { 3, 100 },
            { 4, 120 }
        };
        // Точка, до которой считаем по схеме Горнера
        bmp_real x_div = bmp_real(-0.1);
        BmpVector I;
        for (int i = 0; i < x.size(); i++) {
            auto value = x[i] > x_div ? fdsf::richardson_method(x[i], K_T_MAP[k], k)
                                         : fdsf::Gorner(x[i], k);
                //I.push_back(GornerSchemeForPrecesionY( x0[i], N));
            I.push_back(value);
            std::cout << "x0: " << x[i] << " I_base: " << I[i] << std::endl;
        }
        return I;
    }

    void GetValue_w(BmpVector &I_base,
                    BmpVector &I_additional,
                    BmpVector y0,
                    BmpVector x0,
                    BmpVector Y,
                    BmpVector X, bmp_real k) {
        bmp_real y;
        for (size_t i = 0; i < I_base.size(); i++) {
            y = log(1 + exp(x0[i]));
            I_base[i] = pow((I_base[i] * exp(x0[i]) / y0[i]), 1.0 / k);
        }

        for (size_t i = 0; i < I_additional.size(); i++) {
            y = log(1 + exp(X[i]));
            I_additional[i] = pow((I_additional[i] * exp(X[i]) / Y[i]), 1 / k);
        }
    }
}

TEST_CASE("k_1") {
    std::cout.precision(std::numeric_limits<bmp_real>::max_digits10);
    //int N_gorner = 260, k = 1;
    //int N_gorner = 214, k = 2;
    //int N_gorner = 165, k = 3;
    const size_t k = 1;
    // TODO: tests for N = 2 - 6
    const size_t N_base = 2;

    // Расчет значения интеграла в базовых узлах
    BmpVector x0, X, y0, Y;
    fdsf::SetLinearTrigonometricGrid(y0, x0, Y, X, N_base);

    // Расчет схемы Горнера и подсчета интеграла на Гауссовой сетке
    BmpVector I_base = computeIntegral(x0, k);
    BmpVector I_additional = computeIntegral(X, k);
    //GetValue_w(I_base, I_additional, y0, x0, Y, X, k);
}
