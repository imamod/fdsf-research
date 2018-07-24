#include "Common.h"
#include "Fdsf.h"
#include <array>

namespace {

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
    setPreciseOutput();
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
    BmpVector I_base = compute::integer(grid.xByY(y0), k);
    BmpVector I_additional = compute::integer(grid.xByY(Y), k);
    //GetValue_w(I_base, I_additional, y0, x0, Y, X, k);
}
