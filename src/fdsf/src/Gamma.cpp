#include "Gamma.h"
#include "Constants.h"
#include <unordered_map>

// TODO: rename as gamma
BmpReal factorial(BmpReal _k) {
    auto k = static_cast<double>(_k);
    std::unordered_map<double, BmpReal> SUPPORTED_HALFINTEGER_INDICES = {
        { -1.5, -2 * sqrt(pi()) },
        { -0.5, sqrt(pi()) },
        { 0.5, sqrt(pi()) / 2 },
        { 1.5, 3 * sqrt(pi()) / 4 },
        { 2.5, 15 * sqrt(pi()) / 8 },
        { 3.5, 105 * sqrt(pi()) / 16 }
    };

    auto it = SUPPORTED_HALFINTEGER_INDICES.find(k);
    if (it != SUPPORTED_HALFINTEGER_INDICES.end()) {
        return SUPPORTED_HALFINTEGER_INDICES[k];
    }

    return !k ? 1 : BmpReal(k*factorial(k - 1));
}