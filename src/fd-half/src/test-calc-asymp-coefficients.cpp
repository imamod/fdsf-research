#include "Common.h"
#include "ZetaFunction.h"
#include "FileSys.h"

#include <iostream>

namespace {
    // Вычислить коэффициенты асимтотического ряда
    BmpVector calculate(BmpReal k) {
        BmpVector A;
        // TODO: сделать гибко, через limits
        //constexpr int N = 10;
        constexpr int N = 12;

        for (int n = 1; n <= N; ++n) {
            const BmpReal constMember = 2.0 - pow(2, 2 - 2 * n);
            BmpReal prod = 1;
            for (int p = 1; p <= 2 * n; ++p) {
                prod *= k + 2 - p;
            }
            BmpReal result = constMember*zetaFunction(2 * n)*prod;
            A.push_back(result);
        }
        return A;
    }
}

TEST_CASE("coefficients") {
    SECTION("m_3half") {
        BmpVector a = calculate(-1.5);
        json coeff = a;
        filesys::writeFile("assympt_m3half.json", coeff);
    }
    SECTION("m_half") {
        BmpVector a = calculate(-0.5);
        json coeff = a;
        filesys::writeFile("assympt_mhalf.json", coeff);
    }
    SECTION("half") {
        BmpVector a = calculate(0.5);
        json coeff = a;
        filesys::writeFile("assympt_half.json", coeff);
    }
    SECTION("3half") {
        BmpVector a = calculate(1.5);
        json coeff = a;
        filesys::writeFile("assympt_3half.json", coeff);
    }
    SECTION("5half") {
        BmpVector a = calculate(2.5);
        json coeff = a;
        filesys::writeFile("assympt_5half.json", coeff);
    }
    SECTION("7half") {
        BmpVector a = calculate(3.5);
        json coeff = a;
        filesys::writeFile("assympt_7half.json", coeff);
    }
}
