#include "BasicService.h"
#include <functional>

namespace {
    
    struct integration_segment_values {
        size_t n; // Текущий отрезок интегрирования
        size_t N; // Общее число отрезков интегрирования
    };

    using FermiFunction = std::function<BmpReal(const BmpReal& ksi, const BmpReal& x, BmpReal k, BmpReal a, const integration_segment_values& isv)>;

    // for others half-integer k
    BmpReal fermi_dirak_half_integer(BmpReal ksi, BmpReal x, 
                                      BmpReal k, BmpReal a, 
                                      integration_segment_values isv) {
        BmpReal denom = (1 + ksi) * (BmpReal(isv.N - isv.n) / BmpReal(isv.N));
        //BmpReal denom = (1 - ksi * ksi);
        BmpReal exp_ksi = exp(-a * ksi * ksi / denom);

        return (2*pow(a, k + 1) * pow(ksi, 2 * k + 1) * exp_ksi) /
               (pow(denom, k + 2) * (exp_ksi + exp(-x)));
    }

    // for k = -3/2
    BmpReal fermi_dirak_m3half(BmpReal ksi, BmpReal x,
                                BmpReal k, BmpReal a,
                                integration_segment_values isv) {
        BmpReal exp_ksi = exp(-a * ksi * ksi / (1 - ksi * ksi));
        BmpReal sum_exp = exp_ksi + exp(-x);
        BmpReal exp_diff = exp( -a * ksi * ksi / (1 - ksi * ksi) - x);

        return (-4 * exp_diff * sqrt(abs(a)) * exp_ksi) /
            (pow(1 - ksi * ksi, BmpReal( 3.0 / 2 )) * sum_exp * sum_exp);
    }

    // Euler-Macloren Formulas
    BmpReal trapz(FermiFunction f, BmpReal x, const BmpReal k, size_t N, BmpReal a) {
        BmpReal h = BmpReal(1.0 / N);
        integration_segment_values isv = { 0, N };
        BmpReal u0 = f(0, x, k, a, isv);
        // uN принудительно задаем нулем, чтобы не было переполнения
        BmpReal I = u0 / 2;
        // true work
        for (size_t i = 1; i < N; i++) {
            isv.n = i;
            I += f(i * h, x, k, a, isv);
        }
#if 0
        // TODO: расчет на 2 узлах одновременно??? проверить
        for (size_t i = 1; i < N / 2; i = i + 2) {
            I += f(i * h, x, k, a) + f((N - i) * h, x, k, a);
        }

        if (N == 2) {
            I += f(N * h / 2, x, k, a);
        }
#endif
        return h * I;
    }
}

namespace fdsf {

    BmpReal euler_maclaurin_method(BmpReal x, const BmpReal k, int N) {
        BmpReal a = NewtonsMethod(x, k);
        FermiFunction f = (k == -3.0 / 2) ? fermi_dirak_m3half :
                          fermi_dirak_half_integer;
        return trapz(f, x, k, N, a);
    }

}; //fdsf
