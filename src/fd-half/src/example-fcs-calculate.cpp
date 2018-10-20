#include "Common.h"
#include "FullyConvergedSeries.h"

TEST_CASE("calculate") {
    setPreciseOutput();
    SECTION("m12") {
        INFO("���������� �������� ������� �� ������� k = -1/2 ��� x <=0");
        BmpReal I_0 = fcs::calculate(-0.5, 0);
        std::cout << "k = -0.5 I(0) = " << I_0 << std::endl;
    }
    SECTION("12") {
        INFO("���������� �������� ������� �� ������� k = 1/2 ��� x <=0");
        BmpReal I_0 = fcs::calculate(0.5, 0);
        std::cout << "k = 0.5 I(0) = " << I_0 << std::endl;
    }
    SECTION("32") {
        INFO("���������� �������� ������� �� ������� k = 3/2 ��� x <=0");
        BmpReal I_0 = fcs::calculate(1.5, 0);
        std::cout << "k = 1.5 I(0) = " << I_0 << std::endl;
    }
    SECTION("52") {
        INFO("���������� �������� ������� �� ������� k = 5/2 ��� x <=0");
        BmpReal I_0 = fcs::calculate(2.5, 0);
        std::cout << "k = 2.5 I(0) = " << I_0 << std::endl;
    }
    SECTION("72") {
        INFO("���������� �������� ������� �� ������� k = 7/2 ��� x <=0");
        BmpReal I_0 = fcs::calculate(3.5, 0);
        std::cout << "k = 3.5 I(0) = " << I_0 << std::endl;
    }
}