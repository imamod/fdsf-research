#include <FDSFHalf.h>
#include <boost\multiprecision\cpp_int.hpp>

namespace fdsf {

    bmp_real FFD_half(bmp_real ksi, bmp_real x, bmp_real k)
    {
        // для полуцелых
        // return 2 * pow(t, 2 * k + 1) / (1 + exp(t*t - x));
       // bmp_real ksi = t; // TODO: ksi as param, 0 <= ksi <=1
        return (2*pow(ksi, 2*k+1)*exp(-ksi*ksi/(1-ksi*ksi))) / (pow(1-ksi*ksi, 2*k+2)*(exp(-ksi*ksi/(1-ksi*ksi))+exp(-x)));
    }

    // Euler-Macloren Formulas

    static bmp_real trapz(bmp_real(*Ft)(bmp_real, bmp_real, bmp_real), 
                          bmp_real x, const bmp_real k, int N)
    {

        bmp_real I = 0;
        bmp_real h = bmp_real(1.0/N);
        bmp_real u0 = Ft(0, x, k);
       // bmp_real uN = Ft(1, x, k);
        bmp_real uN = 0;
        I += (u0 + uN) / 2;
        for (size_t i = 1; i < N; i++) {
            I += Ft(i*h, x, k);
        }

        I *= h;

        return I;
    }

    static bmp_real quad(bmp_real(*Ft)(bmp_real, bmp_real, bmp_real),
                         bmp_real x, const bmp_real k, int N)
    {

        bmp_real I = 0;
        bmp_real h = bmp_real(1.0 / N);
        for (size_t i = 0; i < N; i++) {
            I += Ft((i+1/2)*h, x, k);
        }
        I *= h;

        return I;
    }

    bmp_real EM_Simpson(bmp_real x, const bmp_real k, int N)
    {
        return (2*quad(FFD_half, x, k, N) + trapz(FFD_half, x, k, N))/3;
    }

}; //fdsf