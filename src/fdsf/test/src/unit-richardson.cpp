#include "Common.h"
#include "Fdsf.h"
#include "Gamma.h"
#include "Richardson.h"
#include "Print.h"

namespace {

    using FermiFunction = std::function<double(double x, double k, double tau)>;

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

    // Euler-Macloren Formulas
    double trapz(FermiFunction f, double x, const double k, size_t N) {
        double h = 12.0 / N;
        std::vector<double> u;
        for (size_t i = 0; i < N + 1; ++i) {
            u.push_back(f(i * h, x, k));
        }
        u.at(0) /= 2.0;
        u.at(N) /= 2.0;
        double I = 0.0;
        for (auto it = u.rbegin(); it != u.rend(); ++it) {
            I += *it;
        }
        //double I = std::accumulate(u.rbegin(), u.rend(), 0.0);
        return h*I;
    }

    void check(double expected, std::function<double(double, size_t)> f, size_t initialGrid) {
        Richardson r(0, f, initialGrid);
        CHECK_NOTHROW(r.calculate());
        CHECK(abs(expected - 2*r.get()) < 1e-15);
    }

    double f_km12(double x, size_t N) {
        return trapz(fd_m12, x, fdsf::index::M1_HALF, N);
    }

    double f_k12(double x, size_t N) {
        return trapz(fd_m12, x, fdsf::index::P1_HALF, N);
    }

    double f_k32(double x, size_t N) {
        return trapz(fd_m12, x, fdsf::index::P3_HALF, N);
    }

    double f_k52(double x, size_t N) {
        return trapz(fd_m12, x, fdsf::index::P5_HALF, N);
    }

    double f_k72(double x, size_t N) {
        return trapz(fd_m12, x, fdsf::index::P7_HALF, N);
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
        const double funcValueAtZero = 0.604898643421630*factorial(-0.5);
        check(funcValueAtZero, f_km12, 12);
    }
    {
        INFO("Проверка функции ФД индекса k = 1/2 в точке x = 0");
        const double funcValueAtZero = 0.765147024625408*factorial(0.5);
        check(funcValueAtZero, f_k12, 12);
    }
    {
        INFO("Проверка функции ФД индекса k = 3/2 в точке x = 0");
        const double funcValueAtZero = 0.867199889012184*factorial(1.5);
        check(funcValueAtZero, f_k32, 12);
    }
    {
        INFO("Проверка функции ФД индекса k = 5/2 в точке x = 0");
        const double funcValueAtZero = 0.927553577773948*factorial(2.5);
        check(funcValueAtZero, f_k52, 12);
    }
    {
        INFO("Проверка функции ФД индекса k = 7/2 в точке x = 0");
        const double funcValueAtZero = 0.961483656632978*factorial(3.5);
        check(funcValueAtZero, f_k72, 12);
    }
}

