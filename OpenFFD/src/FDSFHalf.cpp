#include "FDSFHalf.h"
#include "Newton.h"
#include <iostream>

namespace fdsf {

    typedef bmp_real(*function)(bmp_real ksi, bmp_real x, bmp_real k, bmp_real a);

    bmp_real FFD_half(bmp_real ksi, bmp_real x, bmp_real k, bmp_real a)
    {
        // для полуцелых
        // return 2 * pow(t, 2 * k + 1) / (1 + exp(t*t - x));

        bmp_real exp_ksi = exp(-a*ksi*ksi / (1 - ksi*ksi));

        bmp_real temp = 2 * k + 1;

        return (2*pow(a, k+1)*pow(ksi, 2*k + 1)*exp_ksi) /
               (pow(1 - ksi*ksi, k + 2)*(exp_ksi + exp(-x)));
    }

    // Euler-Macloren Formulas
    static bmp_real trapz(function f, bmp_real x, const bmp_real k, 
                          int N, bmp_real a)
    {
        bmp_real h = bmp_real(1.0/N);
        bmp_real u0 = f(0, x, k, a);
        // uN принудительно задаем нулем, чтобы не было переполнения
        bmp_real I = u0 / 2; 
#if 0
        // true work
        for (size_t i = 1; i < N; i = i + 2) {
            I += f(i*h, x, k, a);
        }
#endif
//#if 0
        // TODO: расчет на 2 узлах одновременно??? проверить
        for (size_t i = 1; i < N / 2; i = i + 2) {
            I += f(i*h, x, k, a) + f((N - i)*h, x, k, a);
        }

        if (N == 2) {
            I += f(N*h/2, x, k, a);
        }
//#endif
        return h*I;
    }

    static bmp_real quad(function f, bmp_real x, const bmp_real k, 
                         int N, bmp_real a)
    {
        bmp_real I = 0;
        bmp_real h = bmp_real(1.0 / N);
        for (size_t i = 0; i < N; i++) {
            I += f((i + 1.0/2)*h, x, k, a);
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
            I += f((i + 1.0 / 2)*h, x, k, a);
        }
        I *= h;

        return I;
    }

    bmp_real EM_Simpson(bmp_real x, const bmp_real k, int N, bmp_real& a)
    {
        //bmp_real a = newton::NewtonsMethod(x, k);
        a = newton::NewtonsMethod(x, k);
        //std::cout << "a = " << a << std::endl;
        //return (2*quad(FFD_half, x, k, N, a) + trapz(FFD_half, x, k, N, a))/3;
        return trapz(FFD_half, x, k, N, a);
        //return quad(FFD_half, x, k, N, a);
    }

}; //fdsf