#include "Common.h"
#include "FileSys.h"

namespace {

    // Расчет коэффициентов при x<0 с помощью ряда по экспонентам
    BmpVector calculate(size_t N) {
        BmpVector a;
        for (size_t n = 2; n < N; ++n) {
            BmpReal value = 0;
            for (size_t p = 1; p < n; ++p) {
                value += pow(p*(n-p), -0.5);
            }
            a.push_back(value/n);
        }
        return a;
    }
}

TEST_CASE("calculate") {
    BmpVector a = calculate(34);
    filesys::writeFile("iffdLeft.txt", a);
}