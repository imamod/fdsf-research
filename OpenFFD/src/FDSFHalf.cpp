#include "FDSFHalf.h"
#include "Newton.h"
#include <iostream>

namespace fdsf {

    typedef bmp_real(*function)(bmp_real ksi, bmp_real x, 
                                bmp_real k, bmp_real a, integration_segment_values isv);


    bmp_real fermi_dirak_half_integer(bmp_real ksi, bmp_real x, 
                                      bmp_real k, bmp_real a, 
                                      integration_segment_values isv)
    {
        bmp_real denom = (1 + ksi) * (bmp_real(isv.N - isv.n) / bmp_real(isv.N));
        //bmp_real denom = (1 - ksi * ksi);
        bmp_real exp_ksi = exp(-a * ksi * ksi / denom);

        return (2*pow(a, k + 1) * pow(ksi, 2 * k + 1) * exp_ksi) /
               (pow(denom, k + 2) * (exp_ksi + exp(-x)));
    }

    bmp_real fermi_dirak_m3half(bmp_real ksi, bmp_real x, bmp_real k, bmp_real a)
    {
        bmp_real exp_ksi = exp(-a * ksi * ksi / (1 - ksi * ksi));
        bmp_real sum_exp = exp_ksi + exp(-x);
        bmp_real exp_diff = exp( -a * ksi * ksi / (1 - ksi * ksi) - x);

        return (-4 * exp_diff * sqrt(a) * exp_ksi) /
            (pow(1 - ksi * ksi, bmp_real( 3.0 / 2 )) * sum_exp * sum_exp);
    }

    // Euler-Macloren Formulas
    static bmp_real trapz(function f, bmp_real x, const bmp_real k, 
                          int N, bmp_real a)
    {
        bmp_real h = bmp_real(1.0 / N);
        integration_segment_values isv = {0, N};
        bmp_real u0 = f(0, x, k, a, isv);
        // uN принудительно задаем нулем, чтобы не было переполнения
        bmp_real I = u0 / 2; 
//#if 0
        // true work
        for (size_t i = 1; i < N; i++) {
            isv.n = i;
            I += f(i * h, x, k, a, isv);
        }
//#endif
#if 0
        // TODO: расчет на 2 узлах одновременно??? проверить
        for (size_t i = 1; i < N / 2; i = i + 2) {
            I += f(i*h, x, k, a) + f((N - i)*h, x, k, a);
        }

        if (N == 2) {
            I += f(N*h/2, x, k, a);
        }
#endif
        return h*I;
    }

    static bmp_real quad(function f, bmp_real x, const bmp_real k, 
                         int N, bmp_real a)
    {
        bmp_real I = 0;
        bmp_real h = bmp_real(1.0 / N);
        for (size_t i = 0; i < N; i++) {
            //I += f((i + 1.0/2)*h, x, k, a);
        }
        I *= h;
        return I;
    }

    static bmp_real simpson(function f, bmp_real x, 
                            const bmp_real k, int N, bmp_real a)
    {
        bmp_real I = 0;
        bmp_real h = bmp_real(1.0 / N);
        for (size_t i = 0; i < N; i++) {
            //I += f((i + 1.0 / 2)*h, x, k, a);
        }
        I *= h;

        return I;
    }

    bmp_real euler_maclaurin_method(bmp_real x, const bmp_real k, int N, bmp_real& a)
    {
        //bmp_real a = newton::NewtonsMethod(x, k);
        a = newton::NewtonsMethod(x, k);
        //std::cout << "a = " << a << std::endl;
        if (k == -3.0 / 2) {
            //return trapz(fermi_dirak_m3half, x, k, N, a);
        } 
        else {
            return trapz(fermi_dirak_half_integer, x, k, N, a);
        }
    }

}; //fdsf