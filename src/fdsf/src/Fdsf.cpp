/*
* Боевая реализация вычисления функций ФД
*/
#include "Fdsf.h"
#include "FdIndex.h"
#include "FullyConvergedSeries.h"
#include "AsymptoticSeries.h"
#include "FdHalfQuadratures.h"
#include "Constants.h"
#include "JsonFields.h"

namespace {
    // Вычисление ФД полуцелых индексов
    BmpReal calculateHalf(BmpReal x, BmpReal k) {
        if (x <= 0) {
            // Всюду сходящийся ряд для x <=0
            return fcs::calculate(k, x);
        } else if (x >= asympt_series::limits(k).x_min) {
            // Асимптотический ряд для x >= x_min
            return asympt_series::calculate(k, x);
        }
        // Квадратуры 0 <= x <= x_min
        nlohmann::json result = quad::calculate(k, x);
        return result[fd::I];
    }

    // Вычисление положительного значени ФД целого индекса (cм. препринт 2, формула 10)
    BmpReal positiveFdInteger(BmpReal k, BmpReal x) {
        BmpReal fd_m = fcs::calculate(k, -x);
        const BmpReal PI = pi();
        if (fdsf::index::ZERO == k) {
            return fd_m + x;
        } else if (fdsf::index::P1 == k) {
            return -fd_m + x*x/2 + PI*PI / 6;
        } else if (fdsf::index::P2 == k) {
            return fd_m + x*x*x/3 + PI*PI * x/ 3;
        } else if (fdsf::index::P3 == k) {
            BmpReal sqr_x = x*x;
            BmpReal sqr_pi = PI*PI;
            return -fd_m + sqr_x*sqr_x / 4 + sqr_pi * sqr_x / 2 + 7 * sqr_pi*sqr_pi / 60;
        } else if (fdsf::index::P4 == k) {
            BmpReal sqr_x = x*x;
            BmpReal sqr_pi = PI*PI;
            return fd_m + x*sqr_x*sqr_x / 5 + 2*sqr_pi * x*sqr_x / 3 + 7 * x * sqr_pi*sqr_pi / 15;
        }
        throw std::invalid_argument("Unsupported k");
    }

    // Вычисление ФД целых индексов
    BmpReal calculateInteger(BmpReal k, BmpReal x) {
        if (x <= 0) {
            // Всюду сходящийся ряд для x <=0
            return fcs::calculate(k, x);
        }
        return positiveFdInteger(k, x);
    }
}

namespace fdsf {

    /* Функции ФД целого индекса */
    BmpReal fd_0(BmpReal x) {
        return calculateInteger(index::ZERO, x);
    }

    BmpReal fd_1(BmpReal x) {
        return calculateInteger(index::P1, x);
    }

    BmpReal fd_2(BmpReal x) {
        return calculateInteger(index::P2, x);
    }

    BmpReal fd_3(BmpReal x) {
        return calculateInteger(index::P3, x);
    }

    BmpReal fd_4(BmpReal x) {
        return calculateInteger(index::P4, x);
    }

    /* Функции ФД полуцелого индекса */
    BmpReal fd_m3half(BmpReal x) {
        return calculateHalf(x, index::M3_HALF);
    }

    BmpReal fd_m1half(BmpReal x) {
        return calculateHalf(x, index::M1_HALF);
    }

    BmpReal fd_1half(BmpReal x) {
        return calculateHalf(x, index::P1_HALF);
    }

    BmpReal fd_3half(BmpReal x) {
        return calculateHalf(x, index::P3_HALF);
    }

    BmpReal fd_5half(BmpReal x) {
        return calculateHalf(x, index::P5_HALF);
    }

    BmpReal fd_7half(BmpReal x) {
        return calculateHalf(x, index::P7_HALF);
    }

    /* Интегральная функция ФД */
    BmpReal fd_J(BmpReal x) {
        if (x <= 0) {
            // Всюду сходящийся ряд для x <=0
            return fcs::calculateJmhalf(x);
        } else if (x >= 50) {
            // Асимптотический ряд для x >= x_min TODO
            throw std::invalid_argument("Not supported");
        }
        // Квадратуры 0 <= x <= x_min
        nlohmann::json result = quad::calculateJmhalf(x);
        return result[fd::I];
    }
}
