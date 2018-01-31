#define CATCH_CONFIG_MAIN
#include "Common.h"

namespace calc_precision {

    BmpReal halfInteger(BmpReal x, BmpReal k) {
        return halfInteger(BmpVector{ x }, k).front();
    }

    BmpVector halfInteger(BmpVector x, BmpReal k) {
        BmpVector I;
        for (size_t i = 0; i < x.size(); ++i) {
            I.push_back(fdsf::richardson_method(x[i], k));
            //std::cout << "x0: " << x[i] << " I: " << I[i] << std::endl;
        }
        return I;
    }
}