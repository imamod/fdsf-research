#include "Richardson.h"

Richardson::Richardson(double x, std::function<double(double, size_t)> f, size_t initialGrid)
    : m_x(x)
    , m_I(0)
    , m_f(f)
    , m_N(initialGrid) {
}

/* Вычислить значение на сгущающихся сетках */
void Richardson::calculate() {
    const double epsilon = 1e-11;
    m_I = m_f(m_x, m_N);
    double stop_criteria = 0;
    do {
        double I_2n = m_f(m_x, 2 * m_N);
        stop_criteria = abs(m_I / I_2n - 1);
        m_I = I_2n;
        m_N = 2 * m_N;
    } while (stop_criteria > epsilon);
}

/* Получить вычисленное значение интеграла */
const double Richardson::get() const {
    return m_I;
}
