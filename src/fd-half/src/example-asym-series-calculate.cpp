#include "Common.h"
#include "FdIndex.h"
#include "AsymptoticSeries.h"
#include "ZetaFunction.h"
#include <iostream>

TEST_CASE("calculate") {
    setPreciseOutput();
    SECTION("m_3half") {
        INFO("Вычисление значения функции ФД индекса k = -3/2 при x >= x_min");
        BmpReal I_x_44 = asympt_series::calculate(fdsf::index::M3_HALF, 44);
        std::cout << "k = -1.5 I(44) = " << I_x_44 << std::endl;
        BmpReal I_x_52 = asympt_series::calculate(fdsf::index::M3_HALF, 52);
        std::cout << "k = -1.5 I(52) = " << I_x_52 << std::endl;
    }
    SECTION("m_half") {
        INFO("Вычисление значения функции ФД индекса k = -1/2 при x >= x_min");
        BmpReal I_x_min = asympt_series::calculate(fdsf::index::M1_HALF, 39);
        std::cout << "k = -0.5 I(39) = " << I_x_min << std::endl;
    }
    SECTION("half") {
        INFO("Вычисление значения функции ФД индекса k = 1/2 при x >= x_min");
        BmpReal I_x_min = asympt_series::calculate(fdsf::index::P1_HALF, 35);
        std::cout << "k = 0.5 I(35) = " << I_x_min << std::endl;
    }
    SECTION("3half") {
        INFO("Вычисление значения функции ФД индекса k = 3/2 при x >= x_min");
        BmpReal I_x_min = asympt_series::calculate(fdsf::index::P3_HALF, 33);
        std::cout << "k = 1.5 I(33) = " << I_x_min << std::endl;
    }
    SECTION("5half") {
        INFO("Вычисление значения функции ФД индекса k = 5/2 при x >= x_min");
        BmpReal I_x_min = asympt_series::calculate(fdsf::index::P5_HALF, 30);
        std::cout << "k = 2.5 I(30) = " << I_x_min << std::endl;
    }
    SECTION("7half") {
        INFO("Вычисление значения функции ФД индекса k = 7/2 при x >= x_min");
        BmpReal I_x_min = asympt_series::calculate(fdsf::index::P7_HALF, 29);
        std::cout << "k = 3.5 I(29) = " << I_x_min << std::endl;
    }
}

TEST_CASE("zeta") {
    for (size_t n = 2; n <= 24; n += 2) {
        std::cout << n / 2 << " "<< zetaFunction(n) << std::endl;
    }
}

namespace {
    /*
    * Проверка точности asymptotic и квадратурами в точках 30, 35, 40, 45, 50, 55, 60
    */
    void checkAccuracy(BmpReal k, const std::vector<BmpReal>& expected) {
        std::cout << "k = " << k << std::endl;
        const std::vector<BmpReal> REPER_DOTS = { 30, 35, 40, 45, 50 };
        for (size_t i = 0; i < REPER_DOTS.size(); ++i) {
            BmpReal I_0 = asympt_series::calculate(k, REPER_DOTS[i]);
            std::cout << "x = " << REPER_DOTS[i] << ": " << abs(I_0 / expected[i] - 1) << std::endl;
        }
    }
}

/*
 * k = -1.5
 *     x = 30: 6.03628e-13
 *     x = 35: 8.43769e-15
 *     x = 40: 6.66134e-15
 *     x = 45: 2.88658e-15
 *     x = 50: 1.11022e-15
 * k = -0.5
 *     x = 30: 9.76996e-15
 *     x = 35: 1.22125e-15
 *     x = 40: 1.11022e-16
 *     x = 45: 0
 *     x = 50: 2.22045e-16
 * k = 0.5
 *     x = 30: 2.22045e-16
 *     x = 35: 4.44089e-16
 *     x = 40: 0
 *     x = 45: 3.33067e-16
 *     x = 50: 0
 * k = 1.5
 *     x = 30: 0
 *     x = 35: 2.22045e-16
 *     x = 40: 2.22045e-16
 *     x = 45: 3.33067e-16
 *     x = 50: 2.22045e-16
 * k = 2.5
 *     x = 30: 0
 *     x = 35: 2.22045e-16
 *     x = 40: 8.88178e-16
 *     x = 45: 7.77156e-16
 *     x = 50: 4.44089e-16
 * k = 3.5
 *     x = 30: 1.11022e-16
 *     x = 35: 1.11022e-16
 *     x = 40: 2.22045e-16
 *     x = 45: 2.22045e-16
 *     x = 50: 6.66134e-16
 */


TEST_CASE("accuracy") {
    SECTION("m3half") {
        const std::vector<BmpReal> EXPECTED = {
            -0.3656546825930124,
            -0.33840502628678043,
            -0.31647315842739643,
            -0.2983249512504522,
            -0.28298285818087593,
            //-0.2697902993693615,
            //-0.2582876225442183,
        };
        checkAccuracy(fdsf::index::M3_HALF, EXPECTED);
    }
    SECTION("m1half") {
        const std::vector<BmpReal> EXPECTED = {
            10.949421304406611,
            11.828173306266192,
            12.645850688497948,
            13.413677426392377,
            14.139805291011026,
            //14.830377690490632,
            //15.490161585182172,
        };
        checkAccuracy(fdsf::index::M1_HALF, EXPECTED);
    }
    SECTION("half") {
        const std::vector<BmpReal> EXPECTED = {
            109.69481833726648,
            138.1809826590245,
            168.78492259470102,
            201.3687766466349,
            235.81861512588432,
            //272.0382110490754,
            //309.9448732700438,
        };
        checkAccuracy(fdsf::index::P1_HALF, EXPECTED);
    }
    SECTION("3half") {
        const std::vector<BmpReal> EXPECTED = {
            1985.311377746039,
            2913.472994172872,
            4063.3178052469952,
            5450.194657595623,
            7088.5129602482975,
            //8991.89716209947,
            //11173.302913885145,
        };
        checkAccuracy(fdsf::index::P3_HALF, EXPECTED);
    }
    SECTION("5half") {
        const std::vector<BmpReal> EXPECTED = {
            42929.25758509994,
            73324.08962290642,
            116689.92101913081,
            175894.7978015514,
            253992.56857700663,
            //354212.155521165,
            //479948.5008429282,
        };
        checkAccuracy(fdsf::index::P5_HALF, EXPECTED);
    }
    SECTION("7half") {
        const std::vector<BmpReal> EXPECTED = {
            1014417.2017425145,
            2014719.380775691,
            3656385.9176885653,
            6191224.971004079,
            9922878.385070581,
            //15209976.546367215,
            //22469120.83856069,
        };
        checkAccuracy(fdsf::index::P7_HALF, EXPECTED);
    }
}