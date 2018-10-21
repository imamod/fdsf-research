#include "Common.h"
#include "AsymptoticSeries.h"
#include "ZetaFunction.h"
#include <iostream>

TEST_CASE("calculate") {
    setPreciseOutput();
    SECTION("m12") {
        INFO("Вычисление значения функции ФД индекса k = -1/2 при x >= x_min");
        BmpReal I_x_min = asympt_series::calculate(-0.5, 39);
        std::cout << "k = -0.5 I(39) = " << I_x_min << std::endl;
    }
    SECTION("12") {
        INFO("Вычисление значения функции ФД индекса k = 1/2 при x >= x_min");
        BmpReal I_x_min = asympt_series::calculate(0.5, 35);
        std::cout << "k = 0.5 I(35) = " << I_x_min << std::endl;
    }
    SECTION("32") {
        INFO("Вычисление значения функции ФД индекса k = 3/2 при x >= x_min");
        BmpReal I_x_min = asympt_series::calculate(1.5, 33);
        std::cout << "k = 1.5 I(33) = " << I_x_min << std::endl;
    }
    SECTION("52") {
        INFO("Вычисление значения функции ФД индекса k = 5/2 при x >= x_min");
        BmpReal I_x_min = asympt_series::calculate(2.5, 30);
        std::cout << "k = 2.5 I(30) = " << I_x_min << std::endl;
    }
    SECTION("72") {
        INFO("Вычисление значения функции ФД индекса k = 7/2 при x >= x_min");
        BmpReal I_x_min = asympt_series::calculate(3.5, 29);
        std::cout << "k = 3.5 I(29) = " << I_x_min << std::endl;
    }
}

TEST_CASE("zeta") {
    for (size_t n = 2; n <= 24; n += 2) {
        std::cout << n / 2 << " "<< zetaFunction(n) << std::endl;
    }
}