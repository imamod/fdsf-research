#include "Common.h"
#include "Trapz.h"
#include "Gamma.h"

namespace {

    const double epsilon = 1e-11;

    // Подынтегральная функция для индекса k > -3/2
    double fd_m12(double tau, double x, double k) {
        double denom = 1 + exp(tau * tau - x);
        return pow(tau, 2 * k + 1) / denom;
    }

    void check(double fdIndex, double expected) {
        const int N_init = 12;
        int N = N_init;
        double stop_criteria;
        FermiDirakFunction fd = { fdIndex, 0, fd_m12 };
        double I_n = TrapzFD(fd, N_init).trapz(0);
        do {
            double I_2n = TrapzFD(fd, 2*N).trapz(I_n);
            stop_criteria = abs(I_n / I_2n - 1);
            I_n = I_2n;
            N = 2 * N;
        } while (stop_criteria > epsilon);
        const double normalizedValue = 2 * I_n / factorial(fdIndex);
        CHECK(abs(expected - normalizedValue) < 1e-15);
    }
}

TEST_CASE("check") {
    {
        INFO("Вычисление методом трапеций индекса k = -1/2 в точке x = 0");
        check(fdsf::index::M1_HALF, etalon_fd_at_zero::M1_HALF);
    }
    {
        INFO("Вычисление методом трапеций индекса k = 1/2 в точке x = 0");
        check(fdsf::index::P1_HALF, etalon_fd_at_zero::P1_HALF);
    }
    {
        INFO("Вычисление методом трапеций индекса k = 3/2 в точке x = 0");
        check(fdsf::index::P3_HALF, etalon_fd_at_zero::P3_HALF);
    }
    {
        INFO("Вычисление методом трапеций индекса k = 5/2 в точке x = 0");
        check(fdsf::index::P5_HALF, etalon_fd_at_zero::P5_HALF);
    }
    {
        INFO("Вычисление методом трапеций индекса k = 7/2 в точке x = 0");
        check(fdsf::index::P7_HALF, etalon_fd_at_zero::P7_HALF);
    }
}

