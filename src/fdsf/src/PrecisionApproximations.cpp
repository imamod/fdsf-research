/**
* Реализация вычисления функций ФД с помощью прецизионных формул
* TODO: Добавить для полуцелых индексов, интегральной функции
*/

#include "PrecisionApproximations.h"
#include "FdIndex.h"
#include "Errors.h"
#include "Gamma.h"

namespace {

    // Возвращает коэффициенты числителя для заданного индекса
    const BmpVector nomParamsByIndex(BmpReal k) {
        if (fdsf::index::P1 == k) {
            return { 1,
                     0.271511313821436278,
                     0.056266123806058763,
                     0.006742074046934569,
                     0.000516950515533321,
                     0.000019477183676577,
            };
        } else if (fdsf::index::P2 == k) {
            return { 1,
                     0.226381636434069856,
                     0.053368433557479886,
                     0.006290475634079521,
                     0.000502322827445298,
                     0.000018937967508806,
            };
        } else if (fdsf::index::P3 == k) {
            return { 1,
                     0.158348214538045596,
                     0.046064514990930811,
                     0.004886137910884147,
                     0.000433673330597152,
                     0.000017343561379589,
            };
        } else if (fdsf::index::P4 == k) {
            return { 1,
                     0.056014879123090215,
                     0.035111795789180087,
                     0.002183438694367233,
                     0.000246486152552295,
                     0.000009222817788667,
            };
        }
        throw UnsupportedFdFunction();
    }

    // Возвращает коэффициенты знаменателя для заданного индекса
    const BmpVector denomParamsByIndex(BmpReal k) {
        if (fdsf::index::P1 == k) {
            return{ 1,
                    0.021511313821435284,
                    0.023110517572972142,
                    0.000366908157736541,
                    0.000061042440873272,
            };
        } else if (fdsf::index::P2 == k) {
            return{ 1,
                0.038881636434069113,
                0.024304399874277445,
                0.000629098532643319,
                0.000065701816194546,
            };
        } else if (fdsf::index::P3 == k) {
            return{ 1,
                    0.012514881204710761,
                    0.026669340700092963,
                    0.000328543109454736,
                    0.000082091078789006,
            };
        } else if (fdsf::index::P4 == k) {
            return{ 1,
                   -0.061172620876911286,
                    0.027996854281614683,
                   -0.000751214829430754,
                    0.000086068074714292,
            };
        }
        throw UnsupportedFdFunction();
    }
}

namespace prec_approx_formula {
    // Вычислить левую аппроксимацю в точке для заданного k
    BmpReal calculate(BmpReal k, BmpReal x) {
        BmpReal y = log(1 + exp(x));
        BmpVector params = nomParamsByIndex(k);
        BmpReal nominator = 0;
        for (size_t i = 0; i < params.size(); ++i) {
            nominator += params.at(i) * pow(y, i);
        }
        params = denomParamsByIndex(k);
        BmpReal denominator = 0;
        for (size_t i = 0; i < params.size(); ++i) {
            denominator += params.at(i) * pow(y, i);
        }
        return factorial(k)*y*pow(nominator / denominator, k);
    }

}
