#include "FdHalfQuadratures.h"
#include "FdIndex.h"
#include "JsonFields.h"
#include "Richardson.h"
#include "Logger.h"

namespace quad {

    // Подынтегральная функция для индекса k = -3/2
    double fd_m3half(double tau, double x, double k) {
        Logger log("fd_m3half");
        double ch_x = cosh((tau * tau - x) / 2);
        log.info("tau = " + std::to_string(tau) + " x = " + std::to_string(x) + " ch(x) = " + std::to_string(ch_x));
        return pow(ch_x, -2);
    }

    // Подынтегральная функция для индекса k > -3/2
    double fd_m12(double tau, double x, double k) {
        double denom = 1 + exp(tau * tau - x);
        return pow(tau, 2 * k + 1) / denom;
    }

    nlohmann::json calculate(double k, double x) {
        const int initialGrid = 12;
        // Для индекса k =-3/2 отличный от других индексов вид подынтегральной функции
        FermiFunction f = (k == fdsf::index::M3_HALF) ? fd_m3half : fd_m12;
        FermiDirakFunction fd = {k, x, f};
        Richardson r(std::make_shared<FermiDirakFunction>(fd), initialGrid);
        r.calculate();
        nlohmann::json object = nlohmann::json::object();
        object[fd::X] = x;
        // Домножаем значение интеграла на коэффициент перед ним ( смотри формулы (30, 34) препринт 2 )
        const BmpReal coeff = (k == fdsf::index::M3_HALF) ? -1 : 2;
        object[fd::I] = coeff * r.get();
        // TODO: fix
        //object[fd::N_MAX] = N / 2;
        //std::cout << object.dump() << std::endl;
        return object;
    }
}