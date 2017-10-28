#include "TestCommon.h"

namespace {
    bmp_real func_cos(bmp_real x) {
        return 1.0 / (2 - cos(x));
    }

    bmp_real func_demo(bmp_real x) {
        size_t n = 2;
        bmp_real a = 1.75;
        return (a*a - 1)*pow(a, n)*cos(n*x) / pow(1 - 2 * a*cos(x) + a*a, 1);
    }

    bmp_real func_exp_sin(bmp_real x) {
        return exp(sin(x));
    }

    bmp_real func_exp_cos(bmp_real x) {
        return exp(cos(x));
    }
}

TEST_CASE("ExpConverge") {
    // Для статьи о сверхстепенной сходимости
    //epc::Richardson(func_demo, 0, fdsf::PI, true);
    epc::Richardson(func_cos, 0, 1);
    epc::Richardson(func_exp_cos, 0, 1);
    epc::Richardson(func_exp_sin, 0, 1);
}
