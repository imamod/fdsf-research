#include "Common.h"

TEST_CASE("multiply") {
    std::cout << "Enter 0 to finish multiple sequence." << std::endl;
    std::cout.precision(std::numeric_limits<BmpReal>::max_digits10 + 4);
    BmpReal prod = 1;
    while (true) {
        BmpReal a = enter::number("multiplier");
        if (!a) {
            break;
        }
        prod *= a;
    }
    std::cout << "Production = " << prod << std::endl;
}
