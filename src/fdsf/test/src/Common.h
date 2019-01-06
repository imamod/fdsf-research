#pragma once

#include "catch.hpp"
#include "BasicService.h"
#include "Constants.h"
#include "Fdsf.h"
#include "Grid.h"
#include "MatrixUtils.h"

const std::vector<BmpReal> HALF_INDICES = { -3.0 / 2, -1.0 / 2, 1.0 / 2, 3.0 / 2, 5.0 / 2, 7.0 / 2 };
const std::vector<BmpReal> INTEGER_INDICES = { 1, 2, 3, 4 };
const std::vector<BmpReal> ALL_INDICES{ -3.0 / 2, -1.0 / 2, 1.0 / 2, 1, 3.0 / 2, 2, 5.0 / 2, 3, 7.0 / 2, 4 };

namespace etalon_fd_at_zero {
    const double M3_HALF = 0.380104812609684;
    const double M1_HALF = 0.604898643421630;
    const double P1_HALF = 0.765147024625408;
    const double P3_HALF = 0.867199889012184;
    const double P5_HALF = 0.927553577773948;
    const double P7_HALF = 0.961483656632978;
}

namespace compute {
    /**
    * Вычислить функцию полуцелого индекса в точке
    */
    BmpReal halfInteger(BmpReal x, BmpReal k);

    /**
    * Вычислить функцию полуцелого индекса на векторе значений
    */
    BmpVector halfInteger(BmpVector x, BmpReal k);

    /**
    * Вычислить функцию целого индекса на векторе значений
    */
    BmpVector integer(const BmpVector& x, size_t k);
}
