#include "BasicService.h"
#include "Constants.h"

#include <iostream>
#include <fstream>

namespace epc {

    // Euler-Macloren Formulas
    static BmpReal trapz(std::function<BmpReal(const BmpReal&)> f, BmpReal a, BmpReal b, size_t N) {
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
        BmpReal a = fdsf::NewtonsMethod(x, k);
        //std::cout << "a = " << a << std::endl;
        BmpReal exp_ksi = exp(-a*ksi*ksi / (1 - ksi*ksi));

        return (2 * pow(a, k + 1)*pow(ksi, 2 * k + 1)*exp_ksi) /
            (pow(1 - ksi*ksi, k + 2)*(exp_ksi + exp(-x)));
    }

    BmpReal Richardson(std::function<BmpReal(const BmpReal&)> f, BmpReal a, BmpReal b, BmpReal _c, BmpReal _p) {
       // BmpReal c = 1.75;
       // BmpReal p = 0;
        BmpReal c = _c;
        BmpReal p = _p;
        std::cout << "distPole = " << log(c) << ", p = " << p << std::endl;
        auto maxModF = [c, p]() {
            BmpReal q = 1;
            return (c + 1) * pow(c, p) / pow(c - 1, 2 * q - 1);
        };
        auto minModF = [c, p]() {
            BmpReal q = 1;
            return (c - 1) * pow(c, p) / pow(c + 1, 2 * q - 1);
        };
        size_t N = 1;
        BmpReal stop_criteria;
        const BmpReal distToPole = log(c);
        BmpReal I;
        do {
            //N = countEvery ? N + 1 : 2 * N;
            I = trapz(f, a, b, N);
            //BmpReal thError = pi()*maxModF() / (exp(N*distToPole) - 1);
            BmpReal thError = pi()*(maxModF() + minModF()) / (2 * (exp(2 * N * distToPole) - 1));
            BmpReal realError = pi() - I;
            stop_criteria = realError;
            BmpReal result = realError / thError;
            std::cout << "N = " << N
                << ": real error / teoretical = " << result
                << ", lg(teoretical) = " << log10(thError) << std::endl;
            N = 2 * N;
        } while (abs(stop_criteria) > 1e-15);
        return I;
    }

} // exponential convergence
