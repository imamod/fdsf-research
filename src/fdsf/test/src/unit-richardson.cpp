#include "Common.h"
#include "FdIndex.h"
#include "Gamma.h"
#include "Richardson.h"
#include "Print.h"

namespace {

    // Подынтегральная функция для индекса k = -3/2
    double fd_m3half(double tau, double x, double k) {
        double ch_x = cosh((tau * tau - x) / 2);
        return pow(ch_x, -2);
    }

    // Подынтегральная функция для индекса k > -3/2
    double fd_m12(double tau, double x, double k) {
        double denom = 1 + exp(tau * tau - x);
        return pow(tau, 2 * k + 1) / denom;
    }

    void check(double expected, double index) {
        const size_t initialGrid = 12;
        FermiDirakFunction fd(0, index);
        Richardson r(initialGrid, fd);
        CHECK(abs(expected - r.calculate().I) < 1e-15);
    }
}

TEST_CASE("fd") {
    // TODO
/*    {
        INFO("Проверка функции ФД индекса k = -3/2 в точке x = 0");
        const double funcValueAtZero = 0.380104812609684*factorial(-1.5);
        check(funcValueAtZero, f_km32, 12);
    }
*/
    {
        INFO("Проверка функции ФД индекса k = -1/2 в точке x = 0");
        const double funcValueAtZero = etalon_fd_at_zero::M1_HALF *factorial(fdsf::index::M1_HALF);
        check(funcValueAtZero, fdsf::index::M1_HALF);
    }
    {
        INFO("Проверка функции ФД индекса k = 1/2 в точке x = 0");
        const double funcValueAtZero = etalon_fd_at_zero::P1_HALF*factorial(fdsf::index::P1_HALF);
        check(funcValueAtZero, fdsf::index::P1_HALF);
    }
    {
        INFO("Проверка функции ФД индекса k = 3/2 в точке x = 0");
        const double funcValueAtZero = etalon_fd_at_zero::P3_HALF*factorial(fdsf::index::P3_HALF);
        check(funcValueAtZero, fdsf::index::P3_HALF);
    }
    {
        INFO("Проверка функции ФД индекса k = 5/2 в точке x = 0");
        const double funcValueAtZero = etalon_fd_at_zero::P5_HALF*factorial(fdsf::index::P5_HALF);
        check(funcValueAtZero, fdsf::index::P5_HALF);
    }
    {
        INFO("Проверка функции ФД индекса k = 7/2 в точке x = 0");
        const double funcValueAtZero = etalon_fd_at_zero::P7_HALF*factorial(fdsf::index::P7_HALF);
        check(funcValueAtZero, fdsf::index::P7_HALF);
    }
}

