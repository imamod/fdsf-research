/*
* Боевая реализация вычисления функций ФД
*/
#include "Fdsf.h"
#include "FdIndex.h"
#include "FullyConvergedSeries.h"
#include "AsymptoticSeries.h"
#include "Errors.h"
#include "FdHalfQuadratures.h"
#include "GlobalApproximations.h"
#include "PrecisionApproximations.h"
#include "Constants.h"
#include "JsonFields.h"

namespace {

    /* Способ вычисления левой части при x < x* */
    enum NegativeCalculateType {
        FCS = 1,
        PREC_APPROX,
    };

    // Вычисление ФД полуцелых индексов по схеме FCS + QUAD + ASYMPT
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

    // Вычисление ФД полуцелых индексов по схеме PREC_LEFT (см. магистерский диплом)
    BmpReal calculateHalfPrecApprox(BmpReal x, BmpReal k) {
        if (x <= 3) {
            // Всюду сходящийся ряд для x <=0
            return prec_approx_formula::calculateLeft(k, x);
        }
        throw FunctionNotImplemented();
    }

    // Вычисление положительного значени ФД целого индекса (cм. препринт 2, формула 10)
    BmpReal positiveFdInteger(BmpReal k, BmpReal x, BmpReal fd_m) {
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

    // Вычисление ФД целых индексов с помощью всюду сходящегося ряда
    BmpReal calculateInteger(BmpReal k, BmpReal x, uint8_t leftCalculationType) {
        BmpReal fd_m;
        bool isNegativeX = x <= 0;
        switch (leftCalculationType) {
            case NegativeCalculateType::FCS: {
                // Всюду сходящийся ряд для x <=0
                fd_m = fcs::calculate(k, isNegativeX ? x : -x);
            }
            case NegativeCalculateType::PREC_APPROX: {
                // Прецизионные аппроксимации
                fd_m = prec_approx_formula::calculateLeft(k, isNegativeX ? x : -x);
            }
            default:
                throw UnknownCalculationMethod();
        }
        if (isNegativeX) {
            return fd_m;
        }
        return positiveFdInteger(k, x, fd_m);
    }
}

namespace fdsf {

    /* Функции ФД целого индекса */
    BmpReal fd_0(BmpReal x) {
        return calculateInteger(index::ZERO, x, NegativeCalculateType::FCS);
    }

    BmpReal fd_1(BmpReal x) {
        return calculateInteger(index::P1, x, NegativeCalculateType::FCS);
    }

    BmpReal fd_2(BmpReal x) {
        return calculateInteger(index::P2, x, NegativeCalculateType::FCS);
    }

    BmpReal fd_3(BmpReal x) {
        return calculateInteger(index::P3, x, NegativeCalculateType::FCS);
    }

    BmpReal fd_4(BmpReal x) {
        return calculateInteger(index::P4, x, NegativeCalculateType::FCS);
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
            throw FunctionNotImplemented();
        }
        // Квадратуры 0 <= x <= x_min
        nlohmann::json result = quad::calculateJmhalf(x);
        return result[fd::I];
    }
}

/* Реализация вычислений функций ФД с помощью прецизионных аппроксимаций */
/* Функции ФД целого индекса */
BmpReal fdsf::prec_approx::fd_0(BmpReal x) {
    throw FunctionNotImplemented();
}

BmpReal fdsf::prec_approx::fd_1(BmpReal x) {
    return calculateInteger(index::P1, x, NegativeCalculateType::PREC_APPROX);
}

BmpReal fdsf::prec_approx::fd_2(BmpReal x) {
    return calculateInteger(index::P2, x, NegativeCalculateType::PREC_APPROX);
}

BmpReal fdsf::prec_approx::fd_3(BmpReal x) {
    return calculateInteger(index::P3, x, NegativeCalculateType::PREC_APPROX);
}

BmpReal fdsf::prec_approx::fd_4(BmpReal x) {
    return calculateInteger(index::P4, x, NegativeCalculateType::PREC_APPROX);
}

/* Функции ФД полуцелого индекса */
BmpReal fdsf::prec_approx::fd_m3half(BmpReal x) {
    return calculateHalfPrecApprox(x, fdsf::index::M3_HALF);
}

BmpReal fdsf::prec_approx::fd_m1half(BmpReal x) {
    return calculateHalfPrecApprox(x, fdsf::index::M1_HALF);
}

BmpReal fdsf::prec_approx::fd_1half(BmpReal x) {
    return calculateHalfPrecApprox(x, fdsf::index::P1_HALF);
}

BmpReal fdsf::prec_approx::fd_3half(BmpReal x) {
    return calculateHalfPrecApprox(x, fdsf::index::P3_HALF);
}

BmpReal fdsf::prec_approx::fd_5half(BmpReal x) {
    return calculateHalfPrecApprox(x, fdsf::index::P5_HALF);
}

BmpReal fdsf::prec_approx::fd_7half(BmpReal x) {
    return calculateHalfPrecApprox(x, fdsf::index::P7_HALF);
}

/* Реализация вычислений функций ФД с помощью глобальных аппроксимаций */
/* Улучшенная асимптотика */

BmpReal fdsf::global_approx::improved_asympt::fd_1(BmpReal x) {
    return global_formula::calculateFive(x, fdsf::index::P1, GlobalFive::ONE_COEFFICIENT);
}

BmpReal fdsf::global_approx::improved_asympt::fd_2(BmpReal x) {
    return global_formula::calculateFive(x, fdsf::index::P2, GlobalFive::ONE_COEFFICIENT);
}

BmpReal fdsf::global_approx::improved_asympt::fd_3(BmpReal x) {
    return global_formula::calculateFive(x, fdsf::index::P3, GlobalFive::ONE_COEFFICIENT);
}

BmpReal fdsf::global_approx::improved_asympt::fd_4(BmpReal x) {
    return global_formula::calculateFive(x, fdsf::index::P4, GlobalFive::ONE_COEFFICIENT);
}

BmpReal fdsf::global_approx::improved_asympt::fd_m3half(BmpReal x) {
    throw FunctionNotImplemented();
}

BmpReal fdsf::global_approx::improved_asympt::fd_m1half(BmpReal x) {
    return global_formula::calculateFive(x, fdsf::index::M1_HALF, GlobalFive::ONE_COEFFICIENT);
}

BmpReal fdsf::global_approx::improved_asympt::fd_1half(BmpReal x) {
    return global_formula::calculateFive(x, fdsf::index::P1_HALF, GlobalFive::ONE_COEFFICIENT);
}

BmpReal fdsf::global_approx::improved_asympt::fd_3half(BmpReal x) {
    return global_formula::calculateFive(x, fdsf::index::P3_HALF, GlobalFive::ONE_COEFFICIENT);
}

BmpReal fdsf::global_approx::improved_asympt::fd_5half(BmpReal x) {
    return global_formula::calculateFive(x, fdsf::index::P5_HALF, GlobalFive::ONE_COEFFICIENT);
}

BmpReal fdsf::global_approx::improved_asympt::fd_7half(BmpReal x) {
    return global_formula::calculateFive(x, fdsf::index::P7_HALF, GlobalFive::ONE_COEFFICIENT);
}

/* Низкотемпературная асимптотика */

BmpReal fdsf::global_approx::low_temp::fd_1(BmpReal x) {
    return global_formula::calculateFive(x, fdsf::index::P1, GlobalFive::LOW_TEMP);
}

BmpReal fdsf::global_approx::low_temp::fd_2(BmpReal x) {
    return global_formula::calculateFive(x, fdsf::index::P2, GlobalFive::LOW_TEMP);
}

BmpReal fdsf::global_approx::low_temp::fd_3(BmpReal x) {
    return global_formula::calculateFive(x, fdsf::index::P3, GlobalFive::LOW_TEMP);
}

BmpReal fdsf::global_approx::low_temp::fd_4(BmpReal x) {
    return global_formula::calculateFive(x, fdsf::index::P4, GlobalFive::LOW_TEMP);
}

BmpReal fdsf::global_approx::low_temp::fd_m3half(BmpReal x) {
    throw FunctionNotImplemented();
}

BmpReal fdsf::global_approx::low_temp::fd_m1half(BmpReal x) {
    return global_formula::calculateFive(x, fdsf::index::M1_HALF, GlobalFive::LOW_TEMP);
}

BmpReal fdsf::global_approx::low_temp::fd_1half(BmpReal x) {
    return global_formula::calculateFive(x, fdsf::index::P1_HALF, GlobalFive::LOW_TEMP);
}

BmpReal fdsf::global_approx::low_temp::fd_3half(BmpReal x) {
    return global_formula::calculateFive(x, fdsf::index::P3_HALF, GlobalFive::LOW_TEMP);
}

BmpReal fdsf::global_approx::low_temp::fd_5half(BmpReal x) {
    return global_formula::calculateFive(x, fdsf::index::P5_HALF, GlobalFive::LOW_TEMP);
}

BmpReal fdsf::global_approx::low_temp::fd_7half(BmpReal x) {
    return global_formula::calculateFive(x, fdsf::index::P7_HALF, GlobalFive::LOW_TEMP);
}

/* Наилучшая точность */

BmpReal fdsf::global_approx::best_prec::fd_1(BmpReal x) {
    return global_formula::calculateFive(x, fdsf::index::P1, GlobalFive::BEST_PREC);
}

BmpReal fdsf::global_approx::best_prec::fd_2(BmpReal x) {
    return global_formula::calculateFive(x, fdsf::index::P2, GlobalFive::BEST_PREC);
}

BmpReal fdsf::global_approx::best_prec::fd_3(BmpReal x) {
    return global_formula::calculateFive(x, fdsf::index::P3, GlobalFive::BEST_PREC);
}

BmpReal fdsf::global_approx::best_prec::fd_4(BmpReal x) {
    return global_formula::calculateFive(x, fdsf::index::P4, GlobalFive::BEST_PREC);
}

BmpReal fdsf::global_approx::best_prec::fd_m3half(BmpReal x) {
    throw FunctionNotImplemented();
}

BmpReal fdsf::global_approx::best_prec::fd_m1half(BmpReal x) {
    return global_formula::calculateFive(x, fdsf::index::M1_HALF, GlobalFive::BEST_PREC);
}

BmpReal fdsf::global_approx::best_prec::fd_1half(BmpReal x) {
    return global_formula::calculateFive(x, fdsf::index::P1_HALF, GlobalFive::BEST_PREC);
}

BmpReal fdsf::global_approx::best_prec::fd_3half(BmpReal x) {
    return global_formula::calculateFive(x, fdsf::index::P3_HALF, GlobalFive::BEST_PREC);
}

BmpReal fdsf::global_approx::best_prec::fd_5half(BmpReal x) {
    return global_formula::calculateFive(x, fdsf::index::P5_HALF, GlobalFive::BEST_PREC);
}

BmpReal fdsf::global_approx::best_prec::fd_7half(BmpReal x) {
    return global_formula::calculateFive(x, fdsf::index::P7_HALF, GlobalFive::BEST_PREC);
}
