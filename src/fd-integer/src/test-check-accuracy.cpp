#include "Common.h"
#include <iostream>

namespace {

    const BmpReal X_LEFT = 0;

    void checkAccuracy(BmpReal k, BmpReal x, BmpReal left, BmpReal right) {
        setPreciseOutput();
        std::cout << "k = " << k << ", x = " << x << ", d = " << abs(left/right - 1) << std::endl;
    }
}

TEST_CASE("check_integer_index") {
    SECTION("0") {
        // TODO: fix
        BmpReal I_fcs;
        BmpReal I_prec = log(2);
        //checkAccuracy(0, X_LEFT, I_fcs, I_prec);
    }
    SECTION("1") {
        BmpReal I_fcs = 0.82246703342411309;
        BmpReal I_prec = pi()*pi() / 12;
        checkAccuracy(1, X_LEFT, I_fcs, I_prec);
    }
    SECTION("3") {
        BmpReal I_fcs = 5.6821969769834748;
        BmpReal I_prec = 7 * pow(pi(), 4) / 120;
        checkAccuracy(3, X_LEFT, I_fcs, I_prec);
    }
}
