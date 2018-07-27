#include "BasicService.h"

namespace fdsf {
    using Function = std::function<BmpReal(BmpReal x, BmpReal a, BmpReal k)>;

    // Целевая функция
    static BmpReal F(BmpReal x, BmpReal a, BmpReal k) {
        //return a - 3 * (x + log((8 * a + 9.0) / (8 * a - 9.0)));
        BmpReal f = (k == -3.0 / 2) ? BmpReal(a - (x + log((a + 3.0 / 8) / (a - 3.0 / 8))))
            : BmpReal(a - 3 * (k + 7.0 / 8) * (1 + exp(x - a / 3)));
        return f;
    }

    //Первая производная
    static BmpReal FirstDerivative(BmpReal x, BmpReal a, BmpReal k) {
        //return 1 + 432.0 / (64 * a * a - 81);
        BmpReal df = (k == -3.0 / 2) ? BmpReal(1 + 3.0 / (4 * a * a - 9.0 / 16))
            : BmpReal(1 + (k + 7.0 / 8) * exp(x - a / 3));
        return df;
    }

    // Метод Ньютона 
    static BmpReal NewtonInternal(Function f, Function df,
                                  BmpReal x, BmpReal a0,
                                  BmpReal k, BmpReal eps) {
        BmpReal a1 = a0 - f(x, a0, k) / df(x, a0, k);
        while (abs(a0 - a1) > eps) {
            a0 = a1;
            a1 = a0 - f(x, a0, k) / df(x, a0, k);
        }
        return a1;
    }

    BmpReal NewtonsMethod(BmpReal x, BmpReal k) {
        BmpReal a0 = (k == -3.0 / 2) ? 1.13 : BmpReal(3*(k + 7.0 / 8)); // Начальное приближение
        return NewtonInternal(F, FirstDerivative, x, a0, k, 0.001);
    }

} // namespace fdsf
