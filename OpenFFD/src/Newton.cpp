#include "Newton.h"
//#include <cmath.h>

namespace newton {

        typedef bmp_real(*function)(bmp_real x, bmp_real a, bmp_real k);

        // Целевая функция
        static bmp_real Function(bmp_real x, bmp_real a, bmp_real k)
        {
            return a - 3 * (k + 7.0 / 8) * (1 + exp(x - a / 3));
        }

        //Первая производная
        static bmp_real FirstDerivative(bmp_real x, bmp_real a, bmp_real k)
        {
            return 1 + (k + 7.0 / 8) * exp(x - a / 3);
        }

        // Метод Ньютона 
        static bmp_real NewtonInternal(function f, function df,
            bmp_real x, bmp_real a0,
            bmp_real k, bmp_real eps)
        {
            bmp_real a1 = a0 - f(x, a0, k) / df(x, a0, k);
            while (abs(a0 - a1) > eps) {
                a0 = a1;
                a1 = a0 - f(x, a0, k) / df(x, a0, k);
            }

            return a1;
        }

        bmp_real NewtonsMethod(bmp_real x, bmp_real k)
        {
            bmp_real a0 = 3 * (k + 7.0 / 8); // Начальное приближение
            return NewtonInternal(Function, FirstDerivative, x, a0, k, 0.001);
        }

} // newton
