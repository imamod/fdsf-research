#include "FdHalfQuadratures.h"
#include "FdIndex.h"
#include "JsonFields.h"
#include "Richardson.h"
#include "Logger.h"

namespace quad {

    // Вычислить значение ФД индека k для 0 <= x <= x_min
    nlohmann::json calculate(double k, double x) {
        const int initialGrid = 12;
        FermiDirakFunction fd(x, k);
        RichardsonResult result = Richardson(initialGrid, fd).calculate();
        nlohmann::json object = nlohmann::json::object();
        object[fd::X] = x;
        object[fd::I] = result.I;
        // TODO: fix
        object[fd::N_MAX] = result.N / 2;
        //std::cout << object.dump() << std::endl;
        return object;
    }

    // Вычислить значение интегральной ФД для 0 <= x <= x_min
    BmpReal calculateJmhalf(BmpReal x) {
        const int initialGrid = 12;
        FermiDirakFunction fd(fdsf::index::JM2_HALF, exp(x));
        Richardson r(initialGrid, fd);
        return r.calculate().I;
    }
}