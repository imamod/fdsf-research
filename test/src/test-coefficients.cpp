#include "Common.h"
#include <iostream>

/**
 * всюду сходящийся ряд, расчет b_n_k
 * Формулы Эйлера-Маклорена(формула трапеций), равномерная сетка
 */
namespace {
    BmpReal f(size_t n, const BmpReal& tau) {
        BmpReal mul = pow(1 - 2*exp(-tau*tau), n);
        return mul*exp(-tau*tau);
    }

    BmpReal trapz(size_t n, const BmpReal& P) {
        BmpReal F = 0;
        BmpReal h = 8.0 / P;
        for (size_t p = P - 1; p > 0; --p) {
            F += f(n, 8*p/P);
        }
        return (2*F + pow(-1, n))*h;
    }
}

TEST_CASE("coef") {
    const size_t N = 1;
    BmpReal P0 = 8.0;
    BmpReal I2n, In, stopCriteria;
    std::cout.precision(std::numeric_limits<BmpReal>::max_digits10);
    size_t count = 0;
    do {
        In = trapz(N, P0);
        I2n = trapz(N, 2*P0);
        P0 *= 2;
        stopCriteria = abs(I2n - In);
        ++count;
       // std::cout << In << " " << I2n << std::endl;
        std::cout <<count << ": " << stopCriteria << std::endl;
    } while (stopCriteria > 1e-14);
    auto result = I2n / sqrt(fdsf::PI);
    std::cout << "b" << N << " = " << result << std::endl;
}
