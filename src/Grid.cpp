#include "Grid.h"

Grid::Grid(size_t N_base, size_t addNCount)
    : m_N_base(N_base)
    , m_additionalNCount(addNCount)
    , m_base(0)
    , m_additional(0) {}

// Установить сетку базовых узлов по линейному закону
void Grid::setLinearGrid() {
    const BmpReal x_star = BmpReal(3);
    const BmpReal y_star = BmpReal(log(1 + exp(x_star))); // if half-integer
    BmpReal baseSize = BmpReal(2 * m_N_base); // if half-integer & fixed a(N+1)
    const BmpReal y_star_inv = 1 / y_star;

    // Задаются базовые узлы интерполяции
    for (size_t j = 1; j <= baseSize; j++) {
        m_base.push_back(j*y_star_inv / baseSize);
    }

    // Задаются дополнительные точки
    m_additional.push_back(m_base.front() / m_additionalNCount);

    for (size_t i = 1; i < m_additionalNCount; i++) {
        m_additional.push_back(m_additional[i - 1] + m_base.front() / m_additionalNCount);
    }

    for (size_t index = 1; index < m_base.size(); index++) {
        for (size_t i = 0; i < m_additionalNCount; i++) {
            m_additional.push_back(m_additional.back() + (m_base[index] - m_base[index - 1]) / m_additionalNCount);
        }
    }

    // Разворачиваем y
    std::reverse(m_base.begin(), m_base.end());
    for (size_t j = 0; j < baseSize; j++) {
        m_base[j] = 1.0 / m_base[j];
    }

    std::reverse(m_additional.begin(), m_additional.end());
    //std::cout << Y.size() << std::endl;
    for (size_t j = 0; j < m_additional.size(); j++) {
        m_additional[j] = 1.0 / m_additional[j];
    }
}

// Установить сетку базовых узлов по линейно-тригонометрическому закону
void Grid::setLinearTrigonometricGrid() {
    const BmpReal alpha = 2 / (2 + fdsf::PI);
    const BmpReal one = BmpReal(1);
    const BmpReal num2 = BmpReal(2); //if integer
    const BmpReal x_star = BmpReal(3);
    const BmpReal y_star = BmpReal(log(1 + exp(x_star))); // if half-integer
    BmpReal baseSize = BmpReal(2 * m_N_base + 1); // if integer || half-integer & !fixed a(N+1)
    //BmpReal baseSize = BmpReal(2 * N_base ); // if half-integer & fixed a(N+1)

    // Задаются базовые узлы интерполяции
    for (size_t j = 1; j <= baseSize; j++) {
        m_base.push_back(y_star / num2*(num2 * alpha*j / baseSize
            + (one - alpha)*(one - cos(fdsf::PI*j / baseSize))));
    }

    m_additional.push_back(m_base.front() / m_additionalNCount);

    for (size_t i = 1; i < m_additionalNCount; ++i) {
        m_additional.push_back(m_additional[i - 1] + m_base.front() / m_additionalNCount);
    }

    for (size_t index = 1; index < m_base.size(); ++index) {
        for (size_t i = 0; i < m_additionalNCount; ++i) {
            m_additional.push_back(m_additional.back() + (m_base[index] - m_base[index - 1]) / m_additionalNCount);
        }
    }
}

/**
* Задает правую линейно-тригонометрическую сетку базовых узлов справа, плюс 10(?)
* дополнительных точек между каждой парой базовых узлов. Актуально только для полуцелых индексов
*/
void Grid::setLinearTrigonometricGridRight() {
    m_additional.clear();
    m_base.clear();
    //const BmpReal alpha = 2 / (2 + PI);
    const BmpReal alpha = 0.7;
    const BmpReal one = BmpReal(1);
    const BmpReal num2 = BmpReal(2); //if integer
    const BmpReal x_star = BmpReal(3);
    const BmpReal y_star = BmpReal(log(1 + exp(x_star))); // if half-integer
    //BmpReal baseSize = BmpReal(2 * m_N_base + 1); // if integer || half-integer & !fixed a(N+1)
    BmpReal baseSize = BmpReal(2 * m_N_base); // if half-integer & fixed a(N+1)
    //BmpReal baseSize = BmpReal(m_N_base); // if poly approximation

    //const BmpReal y_star_inv = 1 / (y_star * y_star);
    const BmpReal y_star_inv = 1 / y_star;
    //const BmpReal y_star_inv = 1 / pow(y_star, 0.5);
    //const BmpReal y_star_inv = 1 / pow(y_star, 0.25 );
    //const BmpReal y_star_inv = 1 / pow(y_star, 3.0 / 2);

    // Задаются базовые узлы интерполяции
    for (size_t j = 1; j <= baseSize; j++) {
        m_base.push_back(y_star_inv / num2*(num2 * alpha*j / baseSize
            + (one - alpha)*(one - cos(fdsf::PI*j / baseSize))));
    }

    // Задаются дополнительные точки
    m_additional.push_back(m_base.front() / m_additionalNCount);

    for (size_t i = 1; i < m_additionalNCount; i++) {
        m_additional.push_back(m_additional[i - 1] + m_base.front() / m_additionalNCount);
    }

    for (size_t index = 1; index < m_base.size(); index++) {
        for (size_t i = 0; i < m_additionalNCount; i++) {
            m_additional.push_back(m_additional.back() + (m_base[index] - m_base[index - 1]) / m_additionalNCount);
        }
    }

    // Разворачиваем y
    std::reverse(m_base.begin(), m_base.end());
    for (size_t j = 0; j < baseSize; j++) {
        //y_base[j] = 1.0 / pow(y_base[j], 0.5);
        m_base[j] = 1.0 / m_base[j];
        //y_base[j] = 1.0 / (y_base[j] * y_base[j]);
        //y_base[j] = 1.0 / (pow(y_base[j], 4));
        //y_base[j] = 1.0 / (pow(y_base[j], 2.0 / 3));
    }

    std::reverse(m_additional.begin(), m_additional.end());
    //std::cout << Y.size() << std::endl;
    for (size_t j = 0; j < m_additional.size(); ++j) {
        //Y[j] = 1.0 / pow(Y[j], 0.5);
        m_additional[j] = 1.0 / m_additional[j];
        //Y[j] = 1.0 / (Y[j] * Y[j]);
        //Y[j] = 1.0 / pow(Y[j], 4);
        //Y[j] = 1.0 / pow(Y[j], 2.0 / 3);
    }
}

// Получить массив базовых точек
BmpVector Grid::base() const {
    return m_base;
}

// Получить массив дополнительных точек (вместе с базовыми)
BmpVector Grid::additional() const {
    return m_additional;
}

// Получить x по y
BmpVector Grid::xByY(const BmpVector& y) {
    BmpVector x;
    for (const auto& it : y) {
        x.push_back(log(exp(it) - 1));
    }
    return x;
}
