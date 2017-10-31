#include "BasicService.h"
#include "Fdsf.h"
#include "Newton.h"

#include <iostream>
#include <fstream>

namespace epc {

    // Euler-Macloren Formulas
    static bmp_real trapz(function f, bmp_real a, bmp_real b, size_t N) {
        bmp_real h = bmp_real((b - a) / N);
        bmp_real I = (f(a) + f(b)) / 2;
        //bmp_real I = f(a) / 2;

        for (size_t i = 1; i < N; i++) {
            I += f(a + i*h);
        }

        return h*I;
    }

    static bmp_real func_fermi_dirak_half_integer(bmp_real ksi) {
        bmp_real x = -10.0;
        bmp_real k = 3.0 / 2;
        bmp_real a = fdsf::newton::NewtonsMethod(x, k);
        //std::cout << "a = " << a << std::endl;
        bmp_real exp_ksi = exp(-a*ksi*ksi / (1 - ksi*ksi));

        return (2 * pow(a, k + 1)*pow(ksi, 2 * k + 1)*exp_ksi) /
            (pow(1 - ksi*ksi, k + 2)*(exp_ksi + exp(-x)));
    }

    bmp_real Richardson(function f, bmp_real a, bmp_real b, bool countEvery) {
        size_t N = 1;
        bmp_real stop_criteria;
        bmp_real I = trapz(f, a, b, N);
        //bmp_real a_ = 1.75;
        //size_t n = 0;
        //bmp_real I_prec = pow(a_, n)*fdsf::PI / (1 - a_*a_ );
        //bmp_real I_prec = fdsf::PI / ((a_*a_ - 1)*pow(a_, n));
        //bmp_real I_prec = fdsf::PI;
        std::cout << "N = " << N << ": I = " << I << std::endl;
        std::ofstream fout;
        //fout.open("demo.txt");
        //fout.precision(std::numeric_limits<bmp_real>::max_digits10);
        do {
            bmp_real I_2n;

            if (countEvery) {
                I_2n = trapz(f, a, b, N + 1);
                N = N + 1;
            }
            else {
                I_2n = trapz(f, a, b, 2 * N);
                N = 2 * N;
            }

            stop_criteria = (I / I_2n) - 1;
            //stop_criteria = (I_2n / I_prec) - 1;
            I = I_2n;
            
            // std::cout << "N = " << N << ": I = " << I << std::endl;
            //std::cout << "N = " << N << ": d = " << abs(stop_criteria) << std::endl;
            fout << abs(stop_criteria) << std::endl;
        } while (abs(stop_criteria) > fdsf::epsilon);

        //fout.close();
        return I;
    }

} // exponential convergence
