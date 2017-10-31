#include "Newton.h"

namespace fdsf {
    namespace newton {

        typedef BmpReal(*function)(BmpReal x, BmpReal a, BmpReal k);

        // Целевая функция
        static BmpReal Function(BmpReal x, BmpReal a, BmpReal k)
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
        static BmpReal FirstDerivative(BmpReal x, BmpReal a, BmpReal k)
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
        static BmpReal NewtonInternal(function f, function df,
            BmpReal x, BmpReal a0,
            BmpReal k, BmpReal eps)
        {
            BmpReal a1 = a0 - f(x, a0, k) / df(x, a0, k);
            while (abs(a0 - a1) > eps) {
                a0 = a1;
                a1 = a0 - f(x, a0, k) / df(x, a0, k);
            }

            return a1;
        }

        BmpReal NewtonsMethod(BmpReal x, BmpReal k)
        {
            BmpReal a0;

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