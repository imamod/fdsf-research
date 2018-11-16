#pragma once

#include "catch.hpp"
#include "BasicService.h"
#include "Constants.h"
#include "Fdsf-legacy.h"
#include "Grid.h"
#include "MatrixUtils.h"

#include "Enter.h"
#include "Print.h"

#include <iostream>

const std::vector<BmpReal> HALF_INDICES = { -3.0 / 2, -1.0 / 2, 1.0 / 2, 3.0 / 2, 5.0 / 2, 7.0 / 2 };
const std::vector<BmpReal> INTEGER_INDICES = { 1, 2, 3, 4 };
const std::vector<BmpReal> ALL_INDICES { -3.0 / 2, -1.0 / 2, 1.0 / 2, 1, 3.0 / 2, 2, 5.0 / 2, 3, 7.0 / 2, 4 };

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
