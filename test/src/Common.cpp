#define CATCH_CONFIG_MAIN
#include "Common.h"

namespace compute {
    /**
     * Вычислить функцию полуцелого индекса в точке
     */
    BmpReal halfInteger(BmpReal x, BmpReal k) {
        return halfInteger(BmpVector{ x }, k).front();
    }

    /**
     * Вычислить функцию полуцелого индекса на векторе значений
     */
    BmpVector halfInteger(BmpVector x, BmpReal k) {
        BmpVector I;
        for (size_t i = 0; i < x.size(); ++i) {
            I.push_back(fdsf::richardson_method(x[i], k));
            //std::cout << "x0: " << x[i] << " I: " << I[i] << std::endl;
        }
        return I;
    }

    /**
     * Вычислить функцию целого индекса на векторе значений
     */
    BmpVector integer(const BmpVector& x, size_t k) {
        // std::cout << "      I1(0) : " << fdsf::PI*fdsf::PI / 12 << std::endl;
        /**
        * Соответствие значений t для каждого k, подобрано экспериментально
        * BmpReal t = fdsf::get_T_max(X.at(i), k);
        */
        std::array<BmpReal, 4> T_VALUES = { 60, 75, 100, 120 };
        // Точка, до которой считаем по схеме Горнера
        BmpReal x_div = BmpReal(-0.1);
        BmpVector I;
        for (int i = 0; i < x.size(); i++) {
            auto value = x[i] > x_div ? fdsf::richardson_method(x[i], k, T_VALUES[k - 1])
                : fdsf::Gorner(x[i], k);
            //I.push_back(GornerSchemeForPrecesionY( x0[i], N));
            I.push_back(value);
            std::cout << "x0: " << x[i] << " I_base: " << I[i] << std::endl;
        }
        return I;
    }
}