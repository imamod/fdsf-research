#include "Gamma.h"
#include "Constants.h"
#include "FdIndex.h"

#include <unordered_map>

// TODO: rename as gamma
BmpReal factorial(BmpReal _k) {
    auto k = static_cast<double>(_k);
    std::unordered_map<double, BmpReal> SUPPORTED_HALFINTEGER_INDICES = {
        { fdsf::index::M3_HALF, -2 * sqrt(pi()) },
        { fdsf::index::M1_HALF, sqrt(pi()) },
        { fdsf::index::P1_HALF, sqrt(pi()) / 2 },
        { fdsf::index::P3_HALF, 3 * sqrt(pi()) / 4 },
        { fdsf::index::P5_HALF, 15 * sqrt(pi()) / 8 },
        { fdsf::index::P7_HALF, 105 * sqrt(pi()) / 16 }
    };

    auto it = SUPPORTED_HALFINTEGER_INDICES.find(k);
    if (it != SUPPORTED_HALFINTEGER_INDICES.end()) {
        return SUPPORTED_HALFINTEGER_INDICES[k];
    }

    return !k ? 1 : BmpReal(k*factorial(k - 1));
}