#include "FermiDirakFunction.h"
#include "FdIndex.h"
#include "TrapzFD.h"
#include "TrapzJm2Half.h"
#include "Logger.h"
#include "expTau2.h"

namespace {

    // Подынтегральная функция для индекса k = -3/2
    double fd_m3half(double tau, double x, double k) {
        Logger log("fd_m3half");
        double ch_x = cosh((tau * tau - x) / 2);
        //log.info("tau = " + std::to_string(tau) + " x = " + std::to_string(x) + " ch(x) = " + std::to_string(ch_x));
        return pow(ch_x, -2);
    }

    // Подынтегральная функция для индекса k > -3/2
    double fd_m12(double tau, double x, double k) {
        Logger logger("fd_m1half");
        logger.info("tau = " + std::to_string(tau) + ", x = " + std::to_string(x) + ", k = " + std::to_string(k));
        //double p_tau2 = exp(tau * tau);
        double p_tau2 = EXP_TAU2.at(tau);
        double denom = 1 + p_tau2 * exp(-x);
        return pow(tau, 2 * k + 1) / denom;
    }

    // Подынтегральная функция интегральной ФД
    double fd_Jmhalf(double exp_x, double tau, double theta) {
        // Logger logger("fd_Jmhalf");
        double a = exp(tau*tau);
        double t = exp_x;
        if (tau == theta) {
            return log(1 + t / a) - t / (a + t);
        }
        double b = exp(theta*theta);
        double denom = a - b;
        double num = a*log(1 + t / a) - b*log(1 + t / b);
        return num / denom;
    }
}

FermiDirakFunction::FermiDirakFunction(double x, double index)
    : m_params{ index, x } {}

// Вычислить значение
BmpReal FermiDirakFunction::calculate(BmpReal grid, BmpReal previousGridValue) {
    Logger logger("FermiDirakFunction::calculate");
    if (fdsf::index::JM2_HALF == m_params.index) {
        throw std::exception("not implemented");
        // Интегральная функция
        //std::shared_ptr<Trapz> trapzMethod = std::make_shared<TrapzJm2Half>(fd_Jmhalf);
        // return trapz(trapzMethod, grid, previousGridValue);
    }
    // Для индекса k =-3/2 отличный от других индексов вид подынтегральной функции
    fdsf::SubIntegralFunc func = (fdsf::index::M3_HALF == m_params.index) ? fd_m3half : fd_m12;
    std::shared_ptr<Trapz> trapzMethod = std::make_shared<TrapzFD>(func);
    BmpReal result = trapz(trapzMethod, grid, previousGridValue);
    return result;
}

BmpReal FermiDirakFunction::trapz(const std::shared_ptr<Trapz>& trapz, BmpReal grid, BmpReal previousGridValue) {
    Logger logger("FermiDirakFunction::trapz");
    if (previousGridValue) {
        return trapz->economicGrid(grid, previousGridValue, m_params);
    }
    return trapz->firstGrid(grid, m_params);
}

// Возвращает индекс вычисляемой функции
BmpReal FermiDirakFunction::index() const {
    return m_params.index;
}
