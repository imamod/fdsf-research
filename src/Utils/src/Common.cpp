#define CATCH_CONFIG_MAIN
#include "Common.h"
#include "Fdsf.h"
#include "FdIndex.h"

namespace {

    std::function<BmpReal(BmpReal)> fdFunc(BmpReal k) {
        std::function<BmpReal(BmpReal)> f;
        if (fdsf::index::M3_HALF == k) {
            f = fdsf::fd_m3half;
        } else if (fdsf::index::M1_HALF == k) {
            f = fdsf::fd_m1half;
        } else if (fdsf::index::P1_HALF == k) {
            f = fdsf::fd_1half;
        } else if (fdsf::index::P1 == k) {
            f = fdsf::fd_1;
        } else if (fdsf::index::P3_HALF == k) {
            f = fdsf::fd_3half;
        } else if (fdsf::index::P2 == k) {
            f = fdsf::fd_2;
        } else if (fdsf::index::P5_HALF == k) {
            f = fdsf::fd_5half;
        } else if (fdsf::index::P3 == k) {
            f = fdsf::fd_3;
        } else if (fdsf::index::P7_HALF == k) {
            f = fdsf::fd_7half;
        } else if (fdsf::index::P4 == k) {
            f = fdsf::fd_4;
        }
        return f;
    }
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
            auto f = fdFunc(k);
            I.push_back(f(x[i]));
        }
        return I;
    }

}