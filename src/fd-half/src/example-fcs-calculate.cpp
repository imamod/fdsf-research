#include "Common.h"
#include "FullyConvergedSeries.h"

TEST_CASE("calculate") {
    setPreciseOutput();
    SECTION("m3half") {
        INFO("Вычисление значения функции ФД индекса k = -3/2 при x <=0");
        BmpReal I_0 = fcs::calculate(-1.5, 0);
        std::cout << "k = -1.5 I(0) = " << I_0 << std::endl;
    }
    SECTION("mHalf") {
        INFO("Вычисление значения функции ФД индекса k = -1/2 при x <=0");
        BmpReal I_0 = fcs::calculate(-0.5, 0);
        std::cout << "k = -0.5 I(0) = " << I_0 << std::endl;
    }
    SECTION("0") {
        INFO("Вычисление значения функции ФД индекса k = 0 при x <=0");
        BmpReal I_0 = fcs::calculate(0, 0);
        std::cout << "k = 0 I(0) = " << I_0 << std::endl;
    }
    SECTION("half") {
        INFO("Вычисление значения функции ФД индекса k = 1/2 при x <=0");
        BmpReal I_0 = fcs::calculate(0.5, 0);
        std::cout << "k = 0.5 I(0) = " << I_0 << std::endl;
    }
    SECTION("1") {
        INFO("Вычисление значения функции ФД индекса k = 1 при x <=0");
        BmpReal I_0 = fcs::calculate(1, 0);
        std::cout << "k = 1 I(0) = " << I_0 << std::endl;
    }
    SECTION("3half") {
        INFO("Вычисление значения функции ФД индекса k = 3/2 при x <=0");
        BmpReal I_0 = fcs::calculate(1.5, 0);
        std::cout << "k = 1.5 I(0) = " << I_0 << std::endl;
    }
    SECTION("2") {
        INFO("Вычисление значения функции ФД индекса k = 2 при x <=0");
        BmpReal I_0 = fcs::calculate(2, 0);
        std::cout << "k = 2 I(0) = " << I_0 << std::endl;
    }
    SECTION("5half") {
        INFO("Вычисление значения функции ФД индекса k = 5/2 при x <=0");
        BmpReal I_0 = fcs::calculate(2.5, 0);
        std::cout << "k = 2.5 I(0) = " << I_0 << std::endl;
    }
    SECTION("3") {
        INFO("Вычисление значения функции ФД индекса k = 3 при x <=0");
        BmpReal I_0 = fcs::calculate(3, 0);
        std::cout << "k = 3 I(0) = " << I_0 << std::endl;
    }
    SECTION("7half") {
        INFO("Вычисление значения функции ФД индекса k = 7/2 при x <=0");
        BmpReal I_0 = fcs::calculate(3.5, 0);
        std::cout << "k = 3.5 I(0) = " << I_0 << std::endl;
    }
    SECTION("4") {
        INFO("Вычисление значения функции ФД индекса k = 4 при x <=0");
        BmpReal I_0 = fcs::calculate(4, 0);
        std::cout << "k = 4 I(0) = " << I_0 << std::endl;
    }
}