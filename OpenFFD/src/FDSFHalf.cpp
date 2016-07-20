#include <FDSFHalf.h>
#include <boost\multiprecision\cpp_int.hpp>

namespace fdsf {
    // Метод Ньютона 
    typedef double(*function)(double x, double a, double k);

    double TangentsMethod(function f, function df,
        double x, double an, double k, double eps)
    {
        double a0 = an;
        double a1 = a0 - f(x, a0, k) / df(x, a0, k);
        size_t count = 1;
        while (abs(a0 - a1) > eps) {
            a0 = a1;
            a1 = a0 - f(x, a0, k) / df(x, a0, k);
            //std::cout << count << ": a1 = " << a1 << std::endl;
            count++;
        }

        return a1;
    }

    double MyFunction(double x, double a, double k)
    {
        return a - 3 * (k + 7.0 / 8) * (1 + exp(x - a / 3));
    }

    //Первая производная
    double MyDerivative(double x, double a, double k)
    {
        return 1 + (k + 7.0 / 8) * exp(x - a / 3);
    }

    bmp_real FFD_half(bmp_real ksi, bmp_real x, bmp_real k, bmp_real a)
    {
        // для полуцелых
        // return 2 * pow(t, 2 * k + 1) / (1 + exp(t*t - x));

        bmp_real exp_ksi = exp(-a*ksi*ksi / (1 - ksi*ksi));
        return (2*pow(a, k+1)*pow(ksi, 2*k + 1)*exp_ksi) /
            (pow(1 - ksi*ksi, k + 2)*(exp_ksi + exp(-x)));
    }

    // Euler-Macloren Formulas
    static bmp_real trapz(bmp_real(*Ft)(bmp_real, bmp_real, bmp_real, bmp_real), 
                          bmp_real x, const bmp_real k, int N, bmp_real a)
    {

        bmp_real I = 0;
        bmp_real h = bmp_real(1.0/N);
        bmp_real u0 = Ft(0, x, k, a);
       // bmp_real uN = Ft(1, x, k);
        bmp_real uN = 0;
        I += (u0 + uN) / 2;
        for (size_t i = 1; i < N; i++) {
            I += Ft(i*h, x, k, a);
        }

        I *= h;

        return I;
    }

    static bmp_real quad(bmp_real(*Ft)(bmp_real, bmp_real, bmp_real, bmp_real),
                         bmp_real x, const bmp_real k, int N, bmp_real a)
    {
        bmp_real I = 0;
        bmp_real h = bmp_real(1.0 / N);
        for (size_t i = 0; i < N; i++) {
            I += Ft((i + 1.0/2)*h, x, k, a);
        }
        I *= h;
        return I;
    }

    static bmp_real simpson(bmp_real(*Ft)(bmp_real, bmp_real, bmp_real, bmp_real),
        bmp_real x, const bmp_real k, int N, bmp_real a)
    {
        bmp_real I = 0;
        bmp_real h = bmp_real(1.0 / N);
        for (size_t i = 0; i < N; i++) {
            I += Ft((i + 1.0 / 2)*h, x, k, a);
        }
        I *= h;

        return I;
    }

    bmp_real EM_Simpson(bmp_real x, const bmp_real k, int N)
    {
        //Выбор начального приближения
        bmp_real a0 = 3 * (k + 7.0 / 8);
        //Пример вызова функции
        bmp_real a = TangentsMethod(MyFunction, MyDerivative, x, a0, k, 0.001);
        std::cout << "a = " << a << std::endl;
        //return (2*quad(FFD_half, x, k, N, a) + trapz(FFD_half, x, k, N, a))/3;
        return trapz(FFD_half, x, k, N, a);
        //return quad(FFD_half, x, k, N, a);
    }

}; //fdsf