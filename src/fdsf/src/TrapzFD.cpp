#include "TrapzFD.h"
#include "Logger.h"

TrapzFD::TrapzFD(fdsf::SubIntegralFunc func)
    : Trapz(func) {
    Logger logger("TrapzFD::TrapzFD()");
}

/* Прямая реализация метода трапеци (формула (20) в препринте 2)*/
double TrapzFD::firstGrid(BmpReal initialGrid, const fdsf::Params& params) {
    Logger logger("TrapzFD::firstGrid()");
    BmpReal N = initialGrid;
    fdsf::SubIntegralFunc fd = subIntegralFunction();
    const double h = 12.0 / N;
    // Точка u(N) Берется с весом 1/2
    double I = (fd(N * h, params.x, params.index) / 2);
    // Остальные точки берутся с весом 1
    for (int i = N - 1; i > 0; --i) {
        I += fd(i * h, params.x, params.index);
    }
    // Точка u(0) Берется с весом 1/2
    I += (fd(0, params.x, params.index) / 2);
    return h*I;
}

/* Экономичная реализация метода трапеций (формула (21) в препринте 2) */
double TrapzFD::economicGrid(BmpReal initialGrid, double previousValue, const fdsf::Params& params) {
    Logger logger("TrapzFD::economicGrid()");
    BmpReal N = initialGrid;
    const double h = 12.0 / N; // TODO : h = h_prev/ 2
    double gridSum = 0;
    fdsf::SubIntegralFunc fd = subIntegralFunction();
    for (int i =  N - 1; i > 0; i = i - 2) {
        double value = fd(i * h, params.x, params.index);
        gridSum += value;
    }
    // В препринте 2 опечатка
    return (previousValue / 2) + (h * gridSum);
}
