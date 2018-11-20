#include "FdHalfQuadratures.h"
#include "JsonFields.h"
#include "Logger.h"
// TODO: remove
#include <iostream>
#include <limits>
namespace {
    void print(double x, double I, int N) {
        std::cout.precision(std::numeric_limits<double>::max_digits10);
        std::cout << "x = " << x << " N = " << N << " I = " << I << std::endl;
    }
}

namespace quad {

    using FermiFunction = std::function<double(const double& x, double k, double tau)>;

    // Подынтегральная функция для индекса k = -3/2
    double fd_m3half(double tau, double x, double k) {
        Logger log("fd_m3half");
        double ch_x = cosh((tau * tau - x) / 2);
        log.info("tau = " + std::to_string(tau) + " x = " + std::to_string(x) + " ch(x) = " + std::to_string(ch_x));
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

    double euler_maclaurin(double x, double k, int N) {
        FermiFunction f = (k == -1.5) ? fd_m3half : fd_m12;
        return trapz(f, x, k, N);
    }

    // Критерий останова для формул Эйлера-Маклорена
    const double epsilon = 1e-11;

    nlohmann::json calculate(double k, double x) {
        const int N_init = 12;
        int N = N_init;
        double stop_criteria;
        double I_n = euler_maclaurin(x, k, N);
        print(x, I_n, N);
        do {
            double I_2n = euler_maclaurin(x, k, 2 * N);
            stop_criteria = (I_n / I_2n - 1);
            I_n = I_2n;
            N = 2 * N;
            print(x, I_2n, N);
        } while (abs(stop_criteria) > epsilon);
        nlohmann::json object = nlohmann::json::object();
        object[fd::X] = x;
        // Домножаем значение интеграла на коэффициент перед ним ( смотри формулы (30, 34) препринт 2 )
        BmpReal coeff = (k == -1.5) ? -1 : 2;
        object[fd::I] = coeff * I_n;
        object[fd::N_MAX] = N / 2;
        //std::cout << object.dump() << std::endl;
        return object;
    }
}