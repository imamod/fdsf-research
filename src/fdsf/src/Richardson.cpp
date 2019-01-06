#include "Richardson.h"
#include "Logger.h"

Richardson::Richardson(std::shared_ptr<FermiDirakFunction>& f, size_t initialGrid)
    : m_I(0)
    , m_fd(f)
    , m_N(initialGrid) {
    Logger logger("Richardson::Richardson()");
}

/* Вычислить значение на сгущающихся сетках */
void Richardson::calculate() {
    Logger logger("Richardson::calculate()");
    const double epsilon = 1e-11;
    m_I = TrapzFD(*m_fd.get(), m_N).trapz(0);
    double stop_criteria = 0;
    do {
        double I_2n = TrapzFD(*m_fd.get(), 2 * m_N).trapz(m_I);
        stop_criteria = abs(m_I / I_2n - 1);
        m_I = I_2n;
        m_N = 2 * m_N;
    } while (stop_criteria > epsilon);
}

/* Получить вычисленное значение интеграла */
const double Richardson::get() const {
    return m_I;
}
