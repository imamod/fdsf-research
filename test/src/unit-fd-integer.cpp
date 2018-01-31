#include "Common.h"
#include "Fdsf.h"
#include <array>

namespace {
    // Для целых
    BmpVector computeIntegral(const BmpVector& x, size_t k) {
        // std::cout << "      I1(0) : " << fdsf::PI*fdsf::PI / 12 << std::endl;
        /**
         * Соответствие значений t для каждого k, подобрано экспериментально
         * BmpReal t = fdsf::get_T_max(X.at(i), k);
         */
        std::array<BmpReal, 4> T_VALUES = { 60, 75, 100, 120 };
        // Точка, до которой считаем по схеме Горнера
        BmpReal x_div = BmpReal(-0.1);
        BmpVector I;
        for (int i = 0; i < x.size(); i++) {
            auto value = x[i] > x_div ? fdsf::richardson_method(x[i], k, T_VALUES[k-1])
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
                    BmpVector X, BmpReal k) {
        for (size_t i = 0; i < I_base.size(); i++) {
            I_base[i] = pow((I_base[i] * exp(x0[i]) / y0[i]), 1.0 / k);
        }

        for (size_t i = 0; i < I_additional.size(); i++) {
            I_additional[i] = pow((I_additional[i] * exp(X[i]) / Y[i]), 1 / k);
        }
    }
}

TEST_CASE("k_1") {
    std::cout.precision(std::numeric_limits<BmpReal>::max_digits10);
    //int N_gorner = 260, k = 1;
    //int N_gorner = 214, k = 2;
    //int N_gorner = 165, k = 3;
    const size_t k = 1;
    // TODO: tests for N = 2 - 6
    // Расчет значения интеграла в базовых узлах
    Grid grid(5);
    grid.setLinearTrigonometricGrid();
    BmpVector y0 = grid.base();
    BmpVector Y = grid.additional();

    // Расчет схемы Горнера и подсчета интеграла на Гауссовой сетке
    BmpVector I_base = computeIntegral(grid.xByY(y0), k);
    BmpVector I_additional = computeIntegral(grid.xByY(Y), k);
    //GetValue_w(I_base, I_additional, y0, x0, Y, X, k);
}
