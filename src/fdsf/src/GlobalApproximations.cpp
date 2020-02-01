/**
 * Реализация вычисления функций ФД с помощью глобальных формул
 * TODO: Добавить для индекса -3/2, интегральной функции
*/
#include "GlobalApproximations.h"
#include "FdIndex.h"
#include "Gamma.h"
#include "Constants.h"

namespace {

    BmpVector oneCoefficientParams(BmpReal k) {
        BmpReal c1 = 3*(1 - pow(2, -k)) / k;
        BmpReal c3 = pi()*pi() * (k + 1);
        if (fdsf::index::M1_HALF == k) {
            return{ c1, 5.77, c3 };
        } else if (fdsf::index::P1_HALF == k) {
            return{ c1, 2.4, c3 };
        } else if (fdsf::index::P1 == k) {
            return{ c1, 1.73, c3 };
        } else if (fdsf::index::P3_HALF == k) {
            return{ c1, 1.26, c3 };
        } else if (fdsf::index::P2 == k) {
            return{ c1, 0.98, c3 };
        } else if (fdsf::index::P5_HALF == k) {
            return{ c1, 0.77, c3 };
        } else if (fdsf::index::P3 == k) {
            return{ c1, 0.62, c3 };
        } else if (fdsf::index::P7_HALF == k) {
            return{ c1, 0.51, c3 };
        } else if (fdsf::index::P4 == k) {
            return{ c1, 0.42, c3 };
        }
        throw std::invalid_argument("Unsupported Fermi-Dirak function");
    }

    BmpVector lowTempParams(BmpReal k) {
        BmpReal c3 = pi()*pi() * (k + 1);
        if (fdsf::index::M1_HALF == k) {
            return{ 1.07, 7.27, c3 };
        } else if (fdsf::index::P1_HALF == k) {
            return{ 1.23, 2.83, c3 };
        } else if (fdsf::index::P1 == k) {
            return{ 1.15, 1.99, c3 };
        } else if (fdsf::index::P3_HALF == k) {
            return{ 1.03, 1.47, c3 };
        } else if (fdsf::index::P2 == k) {
            return{ 0.93, 1.11, c3 };
        } else if (fdsf::index::P5_HALF == k) {
            return{ 0.83, 0.86, c3 };
        } else if (fdsf::index::P3 == k) {
            return{ 0.75, 0.69, c3 };
        } else if (fdsf::index::P7_HALF == k) {
            return{ 0.67, 0.56, c3 };
        } else if (fdsf::index::P4 == k) {
            return{ 0.60, 0.47, c3 };
        }
        throw std::invalid_argument("Unsupported Fermi-Dirak function");
    }

    BmpVector bestPrecisionParams(BmpReal k) {
        if (fdsf::index::M1_HALF == k) {
            return{ 1.846, 5.430, 7.166 };
        } else if (fdsf::index::P1_HALF == k) {
            return{ 1.44, 2.47, 16.58 };
        } else if (fdsf::index::P1 == k) {
            return{ 1.28, 1.78, 21.50 };
        } else if (fdsf::index::P3_HALF == k) {
            return{ 1.14, 1.32, 26.60 };
        } else if (fdsf::index::P2 == k) {
            return{ 0.99, 1.02, 31.42 };
        } else if (fdsf::index::P5_HALF == k) {
            return{ 0.87, 0.801, 36.513 };
        } else if (fdsf::index::P3 == k) {
            return{ 0.78, 0.65, 41.45 };
        } else if (fdsf::index::P7_HALF == k) {
            return{ 0.702, 0.528, 46.430 };
        } else if (fdsf::index::P4 == k) {
            return{ 0.63, 0.44, 51.59 };
        }
        throw std::invalid_argument("Unsupported Fermi-Dirak function");
    }

    // Получить коэффициенты в зависимости от аппроксимации
    BmpVector paramsByType(uint8_t approxType, BmpReal k) {
        switch (approxType) {
            case GlobalFive::ONE_COEFFICIENT:
                return oneCoefficientParams(k);
            case GlobalFive::LOW_TEMP:
                return lowTempParams(k);
            case GlobalFive::BEST_PREC:
                return bestPrecisionParams(k);
            default:
                throw std::invalid_argument("Unsupported global 5 param");
        }
    }
}

/**
 * Вычисление функции ФД по глобальным 5-членным формулам
 * @param точка x, в которой нужно вычислить значение функции ФД
 * @param индекс функции ФД
 * @param тип глобальной формулы
 */
BmpReal global_formula::calculateFive(BmpReal x, BmpReal k, uint8_t approxType) {
    BmpReal y = log(1 + exp(x));
    BmpVector c = paramsByType(approxType, k);
    BmpReal polinom = pow(factorial(k + 1), 6 / k)*(1 + c.at(0)*y + c.at(1)*y*y) + c.at(2)*pow(y, 4) + pow(y, 6);
    return y / (k + 1)*pow(polinom, k / 6);
}
