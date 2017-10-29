#include "Newton.h"

namespace fdsf {
    namespace newton {

        typedef bmp_real(*function)(bmp_real x, bmp_real a, bmp_real k);

        // Целевая функция
        static bmp_real Function(bmp_real x, bmp_real a, bmp_real k)
        {
            if ( k == -3.0 / 2 )
            {
                 //return a - 3 * (x + log((8 * a + 9.0) / (8 * a - 9.0)));
                return a - (x + log((a + 3.0 / 8) / (a - 3.0 / 8)));
            }
            else {
                return a - 3 * (k + 7.0 / 8) * (1 + exp(x - a / 3));
            }
        }

        //Первая производная
        static bmp_real FirstDerivative(bmp_real x, bmp_real a, bmp_real k)
        {
            if (k == -3.0 / 2) {
                //return 1 + 432.0 / (64 * a * a - 81);
                return 1 + 3.0 / (4 * a * a - 9.0 / 16);
            }
            else {
                return 1 + (k + 7.0 / 8) * exp(x - a / 3);
            }
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
            bmp_real a0;

            if (k == -3.0 / 2) {
                a0 = 1.13; // Начальное приближение
            }
            else {
                a0 = 3 * (k + 7.0 / 8); // Начальное приближение
            }

            return NewtonInternal(Function, FirstDerivative, x, a0, k, 0.001);
        }

    } // namespace newton
} // namespace fdsf