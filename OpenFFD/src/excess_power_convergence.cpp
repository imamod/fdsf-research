#include "FDSFInteger.h"
//#include "FDSFHalf.h"
#include "excess_power_convergence.h"

#include <iostream>
#include <fstream>

namespace epc {
    typedef bmp_real(*function)(bmp_real x);

    // Euler-Macloren Formulas
    static bmp_real trapz(function f, bmp_real a, bmp_real b, size_t N)
    {
        bmp_real h = bmp_real((b - a) / N);
        bmp_real I = (f(a) + f(b)) / 2;

        for (size_t i = 1; i < N; i++) {
            I += f(a + i*h);
        }

        return h*I;
    }

    static bmp_real func_cos(bmp_real x)
    {
        return 1.0 / (2 - cos(x));
        // return cos(x);
    }

    static bmp_real func_demo(bmp_real x)
    {
        size_t n = 0;
        bmp_real a = 2;
        //return cos(n*x)/(1 - 2*a*cos(x) + a*a);
        return (a*a - 1)*pow(a, n)*cos(n*x) / pow(1 - 2 * a*cos(x) + a*a, 2);
    }

    static bmp_real func_exp_sin(bmp_real x)
    {
        return exp(sin(x));
    }

    static bmp_real func_exp_cos(bmp_real x)
    {
        return exp(cos(x));
    }

    static bmp_real Richardson(function f, bmp_real a, bmp_real b)
    {
        size_t N = 1;
        bmp_real stop_criteria;
        bmp_real I = trapz(f, a, b, N);
        bmp_real a_ = 2;
        size_t n = 0;
        //bmp_real I_prec = pow(a_, n)*fdsf::PI / (1 - a_*a_ );
        //bmp_real I_prec = fdsf::PI / ((a_*a_ - 1)*pow(a_, n));
        //bmp_real I_prec = fdsf::PI;
        std::cout << "N = " << N << ": I = " << I << std::endl;
        std::ofstream fout;
        fout.open("demo.txt");
        fout.precision(std::numeric_limits<bmp_real>::max_digits10);
        do {
            bmp_real I_2n = trapz(f, a, b, N + 1);

            stop_criteria = (I / I_2n) - 1;
            //stop_criteria = (I_2n / I_prec) - 1;
            I = I_2n;
            N = N + 1;
            // std::cout << "N = " << N << ": I = " << I << std::endl;
            std::cout << "N = " << N << ": d = " << abs(stop_criteria) << std::endl;
            fout << abs(stop_criteria) << std::endl;
        } while (abs(stop_criteria) > fdsf::epsilon);

        fout.close();
        return I;
    }

    void checkTrapz(bmp_real a, bmp_real b)
    {
        Richardson(func_demo, a, b);
        //Richardson(func_cos, a, b);
        //Richardson(func_exp_cos, a, b);
    }

} // excess power convergence