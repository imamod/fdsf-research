#include "Common.h"
#include <fstream>

TEST_CASE("calc-half-integer") {
    auto k = enter::index();
    setPreciseOutput();
    std::cout << "Calculate: " << std::endl
              << "¹0 : single value" << std::endl
              << "¹1 : range" << std::endl;
    bool calculateRange = enter::number(2);
    if (calculateRange) {
        BmpReal x_left = enter::number("left bound");
        BmpReal x_right = enter::number("right bound");
        REQUIRE(x_right >= x_left);
        BmpReal span = enter::number("span");
        REQUIRE(span > 0);
        size_t N = floor((x_right - x_left)/span);
        BmpVector x;
        for (size_t i = 0; i < N; ++i) {
            x.push_back(x_left + i*span);
        }
        auto I = compute::halfInteger(x, k);
        std::map<BmpReal, BmpReal> xFx;
        for (size_t i = 0; i < x.size(); ++i) {
            xFx[x[i]] = I[i];
        }
        std::cout << std::endl << "Otput : " << std::endl
                  << "¹0 : console" << std::endl
                  << "¹1 : file" << std::endl;
        bool writeToFile = enter::number(2);
        if (writeToFile) {
            std::ofstream file(enter::string("filename"));
            for (auto const& it : xFx) {
                file.precision(6);
                file << it.first << std::setw(4) << " : ";
                file.precision(std::numeric_limits<BmpReal>::max_digits10);
                file << it.second << std::endl;
            }
            file.close();
        } else {
            print::map(xFx);
        }
    } else {
        BmpReal x = enter::number("x");
        std::cout << compute::halfInteger(x, k) << std::endl;
    }
}
