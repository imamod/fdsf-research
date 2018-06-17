#include "Common.h"
#include <iostream>
#include <map>

namespace {

    const std::map<size_t, BmpReal> BERNULLI = {
        { 0, 1 },
        { 2, 1.0 / 6 },
        { 4, -1.0 / 30 },
        { 6, 1.0 / 42 },
        { 8, -1.0 / 30 },
        { 10, 5.0 / 66 },
        { 12, -691.0 / 2730 },
    };

    /**
    dzeta(2n) = 2^(2n-1)*pi^(2n)*|B_2n| / (2n)!
    */
    
    uint64_t factorial(BmpReal n) {
        uint64_t prod = 1;
        for (size_t i = 1; i <= n; ++i) {
            prod *= i;
        }
        return prod;
    }
    /**
     * return dzeta / PI^(2 * n)
     */
    BmpReal dzetaByBernulli(BmpReal b, BmpReal n) {
        //BmpReal dzeta = 0;
        return pow(2, (2 * n - 1)) * abs(b) / factorial(2 * n);
    }
}

TEST_CASE("dzetaByBernulli") {
    std::cout.precision(std::numeric_limits<BmpReal>::max_digits10);
    for (size_t i = 0; i <= 12; i += 2) {
        std::cout << dzetaByBernulli(BERNULLI.at(i), i) << std::endl;
    }
}
