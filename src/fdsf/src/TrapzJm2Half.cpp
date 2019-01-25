#include "TrapzJm2Half.h"
#include "Logger.h"


TrapzJm2Half::TrapzJm2Half(BmpReal grid, const FermiDirakFunction& fd)
    : Trapz(grid, fd) {}

/*********************************************************
*                       PRIVATE
*********************************************************/

// Веса функции
double TrapzJm2Half::weights(int n, int m) {
    //Logger log("weights");
    //log.info("n = " + std::to_string(n) + " m = " + std::to_string(m));
    // в углу вес g = 1/8
    if ((n == m) && !n) {
        return 1.0 / 8;
    }
    // Половинки на диагонали и вертикали g = 1/2
    if ((n == m) || (!n && m)) {
        return 1.0 / 2;
    }
    // Для остальных вес g = 1
    return 1;
}

/* Прямая реализация метода трапеций (формула (?) в препринте 2)*/
double TrapzJm2Half::firstGrid() {
    // Верхняя граница, пока непонятно, какая будет
    double T = 12;
    BmpReal N = getN();
    FermiDirakFunction fd = getFd();
    const double h = T / N;
    double I = 0.0;
    for (int n = N; n > -1; --n) {
        for (int m = N; m >= n; --m) {
            I += weights(n, m) * fd.func(fd.x, n * h, m * h);
        }
    }
    // Умножаем на h^2, т.к. двойно интеграл
    // Умножаем на 2 , т.к. интегрируем по треугольнику
    return 2 * h*h*I;
}

/* Экономичная реализация метода трапеций (формула (21) в препринте 2) */
double TrapzJm2Half::economicGrid(double previousValue) {
    Logger logger("TrapzJm2Half::economicGrid()");
    BmpReal N = getN();
    FermiDirakFunction fd = getFd();
    const double h = 12.0 / N;
    double gridSum = 0;
    /*
    for (int i = 2 * m_N - 1; i > 0; i = i - 2) {
        double value = m_fd.func(i * h, m_fd.x, m_fd.index);
        gridSum += value;
    }*/
    // В препринте 2 опечатка
    return previousValue / 2 + h * gridSum;
}
