#include "AsymptoticSeries.h"
#include <iostream>

AsymptoticSeries::AsymptoticSeries(BmpReal k, BmpReal x)
    : m_seriesSum(sum(k, x)) {
}

// ПОлучить вычисленное значение суммы для конкретного х
BmpReal AsymptoticSeries::get() const {
    return m_seriesSum;
}

/*******************************************************************************
 *                                  PRIVATE
 ******************************************************************************/

// Рассчитать коэффициенты С
BmpVector AsymptoticSeries::getC(BmpReal k, BmpReal x) {
    BmpVector coefficients;
    BmpReal prod = 1;
    BmpReal kPairsProd = k + 1;
    for (size_t j = 0; j < m_A.size(); ++j) {
        // По асимптотической формуле парное добавление множителей, поэтому далее отнимаем 2
        prod *= kPairsProd*(kPairsProd - 1);
        kPairsProd -= 2;
        coefficients.push_back(2 * (1.0 - pow(2, 1.0 - 2 * (j + 1))) * m_A.at(j) * prod);
    }
    return coefficients;
}

// Вычислить сумму асимптотического ряда
BmpReal AsymptoticSeries::sum(BmpReal k, BmpReal x) {
    BmpVector C = getC(k, x);
    BmpReal seriesSum(1);
    for (size_t j = 0; j < m_A.size(); ++j) {
        seriesSum += pow(x, -2.0*(j + 1))*C.at(j);
    }
    seriesSum *= pow(x, k + 1) / (k + 1);
    return seriesSum;
}
