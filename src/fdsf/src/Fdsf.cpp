/*
* Боевая реализация вычисления функций ФД
*/
#include "Fdsf.h"

#include "FullyConvergedSeries.h"
#include "AsymptoticSeries.h"
#include "FdHalfQuadratures.h"
#include "JsonFields.h"

// TODO: вынести в отдельный модуль
namespace fdsf {
    /* Поддерживаемые индексы функций ФД */
    namespace index {
        const BmpReal M3_HALF = -3.0 / 2;
        const BmpReal M1_HALF = -1.0 / 2;
        const BmpReal P1_HALF = 1.0 / 2;
        const BmpReal P3_HALF = 3.0 / 2;
        const BmpReal P5_HALF = 5.0 / 2;
        const BmpReal P7_HALF = 7.0 / 2;
    }
}

namespace {
    // Вычисление ФД
    BmpReal calculate(BmpReal x, BmpReal k) {
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
}

namespace fdsf {

    /* Функции ФД целого индекса */
     // TODO
    BmpReal fd_0(BmpReal x) {
        throw std::domain_error("Function not implemented");
    }

    BmpReal fd_1(BmpReal x) {
        throw std::domain_error("Function not implemented");
    }

    BmpReal fd_2(BmpReal x) {
        throw std::domain_error("Function not implemented");
    }

    BmpReal fd_3(BmpReal x) {
        throw std::domain_error("Function not implemented");
    }

    BmpReal fd_4(BmpReal x) {
        throw std::domain_error("Function not implemented");
    }

    /* Функции ФД полуцелого индекса */
    BmpReal fd_m3half(BmpReal x) {
        return calculate(x, index::M3_HALF);
    }

    BmpReal fd_m1half(BmpReal x) {
        return calculate(x, index::M1_HALF);
    }

    BmpReal fd_1half(BmpReal x) {
        return calculate(x, index::P1_HALF);
    }

    BmpReal fd_3half(BmpReal x) {
        return calculate(x, index::P3_HALF);
    }

    BmpReal fd_5half(BmpReal x) {
        return calculate(x, index::P5_HALF);
    }

    BmpReal fd_7half(BmpReal x) {
        return calculate(x, index::P7_HALF);
    }
}
