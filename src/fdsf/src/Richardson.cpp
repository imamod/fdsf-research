#include "Richardson.h"
#include "FdIndex.h"
#include "Logger.h"

#include <iostream>

Richardson::Richardson(BmpReal initialGrid, const FermiDirakFunction& fd)
    : m_result{ initialGrid, 0 }
    , m_func(std::make_shared<FermiDirakFunction>(fd)) {
    Logger logger("Richardson::Richardson()");
}

/* Вычислить значение на сгущающихся сетках */
RichardsonResult Richardson::calculate() {
    Logger logger("Richardson::calculate()");
    const double epsilon = 1e-11;
    m_result.I = m_func->calculate(m_result.N);
    // std::cout << "N = " << m_result.N << ", I = " << m_result.I << std::endl;
    double stop_criteria = 0;
    do {
        m_result.N *= 2;
        double I_2n = m_func->calculate(m_result.N, m_result.I);
        stop_criteria = abs(m_result.I / I_2n - 1);
        m_result.I = I_2n;
   //     std::cout << "N = " << m_result.N << ", I = " << m_result.I << std::endl;
    } while (stop_criteria > epsilon);
    // Домножаем значение интеграла на коэффициент перед ним ( смотри формулы (30, 34) препринт 2 )
    const BmpReal coeff = (fdsf::index::M3_HALF == m_func->index()) ? -1 : 2;
    m_result.I *= coeff;
    std::cout << "N = " << m_result.N << ", I = " << m_result.I << std::endl;
    return m_result;
}
