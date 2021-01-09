#include "Common.h"
#include "FdIndex.h"
#include "FullyConvergedSeries.h"
#include <iostream>

TEST_CASE("calculate") {
    setPreciseOutput();
    SECTION("m3half") {
        INFO("Вычисление значения функции ФД индекса k = -3/2 при x <=0");
        BmpReal I_0 = fcs::calculate(fdsf::index::M3_HALF, 0);
        std::cout << "k = -1.5 I(0) = " << I_0 << std::endl;
    }
    SECTION("mHalf") {
        INFO("Вычисление значения функции ФД индекса k = -1/2 при x <=0");
        BmpReal I_0 = fcs::calculate(fdsf::index::M1_HALF, 0);
        std::cout << "k = -0.5 I(0) = " << I_0 << std::endl;
    }
    SECTION("half") {
        INFO("Вычисление значения функции ФД индекса k = 1/2 при x <=0");
        BmpReal I_0 = fcs::calculate(fdsf::index::P1_HALF, 0);
        std::cout << "k = 0.5 I(0) = " << I_0 << std::endl;
    }
    SECTION("3half") {
        INFO("Вычисление значения функции ФД индекса k = 3/2 при x <=0");
        BmpReal I_0 = fcs::calculate(fdsf::index::P3_HALF, 0);
        std::cout << "k = 1.5 I(0) = " << I_0 << std::endl;
    }
    SECTION("5half") {
        INFO("Вычисление значения функции ФД индекса k = 5/2 при x <=0");
        BmpReal I_0 = fcs::calculate(fdsf::index::P5_HALF, 0);
        std::cout << "k = 2.5 I(0) = " << I_0 << std::endl;
    }
    SECTION("7half") {
        INFO("Вычисление значения функции ФД индекса k = 7/2 при x <=0");
        BmpReal I_0 = fcs::calculate(fdsf::index::P7_HALF, 0);
        std::cout << "k = 3.5 I(0) = " << I_0 << std::endl;
    }
}


namespace {
    /*
     * Проверка точности fcs и квадратурами в точках -2, -1, 0
     */
    void checkAccuracy(BmpReal k, const std::vector<BmpReal>& expected) {
        std::cout << "k = " << k << std::endl;
        const std::vector<BmpReal> REPER_DOTS = { -2, -1, 0 };
        for (size_t i = 0; i < REPER_DOTS.size(); ++i) {
            BmpReal I_0 = fcs::calculate(k, REPER_DOTS[i]);
            std::cout << "x = " << REPER_DOTS[i] << ": " << I_0 / expected[i] - 1 << std::endl;
        }
    }
}

/*
 * Проверка точности fcs и квадратурами в точках -2, -1, 0
 * На 03.03.2019 получили результаты
 * k = -1.5
 *     x = -2: 2.22045e-16
 *     x = -1: 2.22045e-16
 *     x = 0: 1.11022e-16
 * k = -0.5
 *     x = -2: 2.22045e-16
 *     x = -1: 2.22045e-16
 *     x = 0: 0
 * k = 0.5
 *     x = -2: 2.22045e-16
 *     x = -1: 2.22045e-16
 *     x = 0: 3.33067e-16
 * k = 1.5
 *     x = -2: 2.22045e-16
 *     x = -1: 1.11022e-16
 *     x = 0: 0
 * k = 2.5
 *     x = -2: 0
 *     x = -1: 2.22045e-16
 *     x = 0: 1.11022e-16
 * k = 3.5
 *     x = -2: 0
 *     x = -1: 0
 *     x = 0: 3.33067e-16
 */
TEST_CASE("accuracy") {
    SECTION("m3half") {
        const std::vector<BmpReal> EXPECTED = {
            -0.40108449328644824,
            -0.8394844152932341,
            -1.3474364777155081,
        };
        checkAccuracy(fdsf::index::M3_HALF, EXPECTED);
    }
    SECTION("m1half") {
        const std::vector<BmpReal> EXPECTED = {
            0.2191916075861797,
            0.5211503831079912,
            1.0721549299401913,
        };
        checkAccuracy(fdsf::index::M1_HALF, EXPECTED);
    }
    SECTION("half") {
        const std::vector<BmpReal> EXPECTED = {
            0.11458782392526307,
            0.2905008961699176,
            0.678093895153101,
        };
        checkAccuracy(fdsf::index::P1_HALF, EXPECTED);
    }
    SECTION("3half") {
        const std::vector<BmpReal> EXPECTED = {
            0.17580098885401288,
            0.4608488062901017,
            1.1528038370883613,
        };
        checkAccuracy(fdsf::index::P3_HALF, EXPECTED);
    }
    SECTION("5half") {
        const std::vector<BmpReal> EXPECTED = {
            0.444554453458763,
            1.1859681754434672,
            3.0825860828374183,
        };
        checkAccuracy(fdsf::index::P5_HALF, EXPECTED);
    }
    SECTION("7half") {
        const std::vector<BmpReal> EXPECTED = {
            1.5649662632360704,
            4.213264071926359,
            11.183716751693321,
        };
        checkAccuracy(fdsf::index::P7_HALF, EXPECTED);
    }
}