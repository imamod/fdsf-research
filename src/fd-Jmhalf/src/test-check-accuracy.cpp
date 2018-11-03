#include "Common.h"

namespace {
    void checkAccuracy(double x, double left, double right) {
        std::cout << "x = " << x << ", d = " << left/right - 1<< std::endl;
    }
}

TEST_CASE("check") {
    SECTION("left") {
        double I_fcs  = 0.78323866983319212;
        double I_quad = 0.783238669833194;
        checkAccuracy(0, I_fcs, I_quad);
    }
}