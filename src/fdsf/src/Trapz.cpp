#include "Trapz.h"
#include "Logger.h"

TrapzFD::TrapzFD(const FermiDirakFunction& fd, double grid)
    : m_fd(fd)
    , m_N(grid) {
    Logger logger("TrapzFD::TrapzFD()");
}

// Вычислить методом трапеций значение интеграла на сетке
double TrapzFD::trapz(double previousGridValue) {
    Logger logger("TrapzFD::trapz()");
    if (previousGridValue == 0) {
        return firstGrid();
    }
    return economicGrid(previousGridValue);
}

/*********************************************************
 *                       PRIVATE
 *********************************************************/

/* Прямая реализация метода трапеци (формула (20) в препринте 2)*/
double TrapzFD::firstGrid() {
    Logger logger("TrapzFD::firstGrid()");
    const double h = 12.0 / m_N;
    // Точка u(N) Берется с весом 1/2
    double I = m_fd.func(m_N * h, m_fd.x, m_fd.index) / 2;
    // Остальные точки берутся с весом 1
    for (int i = m_N-1; i > 0; --i) {
        I += m_fd.func(i * h, m_fd.x, m_fd.index);
    }
    // Точка u(0) Берется с весом 1/2
    I += m_fd.func(0, m_fd.x, m_fd.index) / 2;
    return h*I;
}

/* Экономичная реялизация метода трапеций (формула (21) в препринте 2) */
double TrapzFD::economicGrid(double previousValue) {
    Logger logger("TrapzFD::economicGrid()");
    const double h = 12.0 / m_N;
    double gridSum = 0;
    for (int i = 2 * m_N - 1; i > 0; i = i - 2) {
        double value = m_fd.func(i * h, m_fd.x, m_fd.index);
        gridSum += value;
    }
    // В препринте 2 опечатка
    return previousValue / 2 + h * gridSum;
}
