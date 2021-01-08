#define CATCH_CONFIG_MAIN
#include "Common.h"

#include <array>

namespace fdsf {
    // Сгущение по Ричардсону по сеточно-Гауссову методу
    BmpReal richardson_method(BmpReal x, BmpReal k, BmpReal t=0);
}

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
        }
        return I;
    }

}