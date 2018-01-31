#include "Common.h"

TEST_CASE("calc-half-integer") {
    auto k = enter::index();
    std::cout.precision(std::numeric_limits<BmpReal>::max_digits10);
    std::cout << "Enter mode: " << std::endl
              << "¹0 : Calculate for x" << std::endl
              << "¹1 : Calculate range" << std::endl;
    bool calculateRange = enter::number(2);
    if (calculateRange) {
        BmpReal x_left = enter::number("left bound");
        BmpReal x_right = enter::number("right bound");
        REQUIRE(x_right >= x_left);
        BmpReal span = enter::number("span");
        size_t N = floor((x_right - x_left)/span);
        BmpVector x;
        for (size_t i = 0; i < N; ++i) {
            x.push_back(x_left + i*span);
        }
        auto I = calc_precision::halfInteger(x, k);
        std::map<BmpReal, BmpReal> xFx;
        for (size_t i = 0; i < x.size(); ++i) {
            xFx[x[i]] = I[i];
        }
        print::map(xFx);
    } else {
        BmpReal x = enter::number("x");
        std::cout << calc_precision::halfInteger(x, k) << std::endl;
    }
}
