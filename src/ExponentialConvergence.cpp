#include "BasicService.h"
#include "Fdsf.h"
#include "Newton.h"

#include <iostream>
#include <fstream>

namespace epc {

    // Euler-Macloren Formulas
    static BmpReal trapz(function f, BmpReal a, BmpReal b, size_t N) {
        BmpReal h = BmpReal((b - a) / N);
        BmpReal I = (f(a) + f(b)) / 2;
        //BmpReal I = f(a) / 2;

        for (size_t i = 1; i < N; i++) {
            I += f(a + i*h);
        }

        return h*I;
    }

    static BmpReal func_fermi_dirak_half_integer(BmpReal ksi) {
        BmpReal x = -10.0;
        BmpReal k = 3.0 / 2;
        BmpReal a = fdsf::newton::NewtonsMethod(x, k);
        //std::cout << "a = " << a << std::endl;
        BmpReal exp_ksi = exp(-a*ksi*ksi / (1 - ksi*ksi));

        return (2 * pow(a, k + 1)*pow(ksi, 2 * k + 1)*exp_ksi) /
            (pow(1 - ksi*ksi, k + 2)*(exp_ksi + exp(-x)));
    }

    BmpReal Richardson(function f, BmpReal a, BmpReal b, bool countEvery) {
        size_t N = 1;
        BmpReal stop_criteria;
        BmpReal I = trapz(f, a, b, N);
        //BmpReal a_ = 1.75;
        //size_t n = 0;
        //BmpReal I_prec = pow(a_, n)*fdsf::PI / (1 - a_*a_ );
        //BmpReal I_prec = fdsf::PI / ((a_*a_ - 1)*pow(a_, n));
        //BmpReal I_prec = fdsf::PI;
        std::cout << "N = " << N << ": I = " << I << std::endl;
        std::ofstream fout;
        //fout.open("demo.txt");
        //fout.precision(std::numeric_limits<BmpReal>::max_digits10);
        do {
            BmpReal I_2n;

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
