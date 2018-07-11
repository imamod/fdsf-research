#include "Grid.h"
#include "Constants.h"

Grid::Grid(size_t N_base, size_t addNCount)
    : m_N_base(N_base)
    , m_additionalNCount(addNCount)
    , m_base(0)
    , m_additional(0) {}

/* �������� ������ ������� ����� */
BmpVector Grid::base() const {
    return m_base;
}

/* �������� ������ �������������� ����� (������ � ��������) */
BmpVector Grid::additional() const {
    return m_additional;
}

/* �������� x �� y */
BmpVector Grid::xByY(const BmpVector& y) {
    BmpVector x;
    for (const auto& it : y) {
        x.push_back(log(exp(it) - 1));
    }
    return x;
}

/* ��������� ����� ��������������� ������� */
void Grid::setAdditionalDots() {
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

/* ��������� ������ ���������� y = 1/ksi */
void Grid::changeGrid(BmpVector& v) {
    for (auto& item : v) {
        item = 1 / item;
    }
}

// ���������� ����� ������� ����� �� ��������� ������ (������)
void Grid::setLinearGrid() {
    const BmpReal x_star = BmpReal(3);
    const BmpReal y_star = BmpReal(log(1 + exp(x_star))); // if half-integer
    BmpReal baseSize = BmpReal(2 * m_N_base); // if half-integer & fixed a(N+1)
    const BmpReal y_star_inv = 1 / y_star;

    // �������� ������� ���� ������������
    for (size_t j = 1; j <= baseSize; j++) {
        m_base.push_back(j*y_star_inv / baseSize);
    }

    // �������� �������������� �����
    setAdditionalDots();

    // ������������� y
    std::reverse(m_base.begin(), m_base.end());
    changeGrid(m_base);
    std::reverse(m_additional.begin(), m_additional.end());
    changeGrid(m_additional);
}

// ���������� ����� ������� ����� �� �������-������������������� ������
void Grid::setLinearTrigonometricGrid() {
    const BmpReal alpha = 2 / (2 + pi());
    const BmpReal one = BmpReal(1);
    const BmpReal num2 = BmpReal(2); //if integer
    const BmpReal x_star = BmpReal(5);
    const BmpReal y_star = BmpReal(log(1 + exp(x_star))); // if half-integer
    BmpReal baseSize = BmpReal(2 * m_N_base + 1); // if integer || half-integer & !fixed a(N+1)
    //BmpReal baseSize = BmpReal(2 * N_base ); // if half-integer & fixed a(N+1)

    // �������� ������� ���� ������������
    for (size_t j = 1; j <= baseSize; j++) {
        m_base.push_back(y_star / num2*(num2 * alpha*j / baseSize
            + (one - alpha)*(one - cos(pi()*j / baseSize))));
    }

    // �������� �������������� �����
    setAdditionalDots();
}

/**
* ������ ������ �������-������������������ ����� ������� ����� ������, ���� 10(?)
* �������������� ����� ����� ������ ����� ������� �����. ��������� ������ ��� ��������� ��������
*/
void Grid::setLinearTrigonometricGridRight() {
    m_additional.clear();
    m_base.clear();
    const BmpReal alpha = 2 / (2 + pi());
    //const BmpReal alpha = 0.7;
    const BmpReal one = BmpReal(1);
    const BmpReal num2 = BmpReal(2); //if integer
    //const BmpReal x_star = BmpReal(3);
    //const BmpReal x_star = BmpReal(5);
    const BmpReal x_star = BmpReal(7);
    const BmpReal y_star = BmpReal(log(1 + exp(x_star))); // if half-integer
    //BmpReal baseSize = BmpReal(2 * m_N_base + 1); // if integer || half-integer & !fixed a(N+1)
    BmpReal baseSize = BmpReal(2 * m_N_base); // if half-integer & fixed a(N+1)
    //BmpReal baseSize = BmpReal(m_N_base); // if poly approximation

    //const BmpReal y_star_inv = 1 / (y_star * y_star);
    //const BmpReal y_star_inv = 1 / y_star;
    const BmpReal P(1);
    //const BmpReal P(0.9);
    //const BmpReal P(0.95);
    const BmpReal y_star_new = 1 / pow(y_star, P);

    // �������� ������� ���� ������������
    for (size_t j = 1; j <= baseSize; j++) {
        //m_base.push_back(y_star_inv / num2*(num2 * alpha*j / baseSize
        //    + (one - alpha)*(one - cos(pi()*j / baseSize))));
        /* ������� ���� �������� ������� �� 10.02.2018 */
        /*m_base.push_back(y_star_inv / num2*(num2 * (one - alpha)*j / baseSize
            + alpha*(one - cos(pi()*j / baseSize))));*/
        /* ����� ������� �� �������� p*/
        BmpReal linearPart = num2 * (one - alpha)*j / baseSize;
        BmpReal trigonometricPart = alpha*(one - cos(pi()*j / baseSize));
        m_base.push_back(pow(y_star_new / num2*(linearPart + trigonometricPart), 1 / P));
    }

    // �������� �������������� �����
    setAdditionalDots();

    // ������������� y
    std::reverse(m_base.begin(), m_base.end());
    changeGrid(m_base);
    std::reverse(m_additional.begin(), m_additional.end());
    changeGrid(m_additional);
}

/**
 * ������ ����������. ������������ ����� � �������������� ������������� �����������.
 * ksi = 1/y; delta - ����������, ���������� �� ����� � ksi.
 * eta(i) = ksi(i) + 0.5*(ksi(i+1)-ksi(i-1))*tau*(1+sqrt(delta(i-0.5)/delta(i+0.5)))/(1-sqrt(delta(i-0.5)/delta(i+0.5)))
 * tau = 0.5
 */
void Grid::shiftLinTrigGrid(const BmpVector& delta, BmpReal tau) {
    // ����� �������� �������s �����
    BmpVector ksi;//(m_base);
    for (auto const& it : m_base) {
        ksi.push_back(1 / it);
    }
    std::reverse(ksi.begin(), ksi.end());

    BmpVector eta;
    eta.push_back(ksi.front());
    for (auto i = 1; i < ksi.size() - 1; ++i) {
        auto distance = (ksi[i + 1] - ksi[i - 1]) / 2;
        auto num = 1 - sqrt(delta[i - 1] / delta[i]);
        auto denom = 1 + sqrt(delta[i - 1] / delta[i]);
        eta.push_back(ksi[i] + tau*distance*num / denom);
    }
    eta.push_back(ksi.back());

    // Magic for correct filling
    m_base.clear();
    m_base = eta;

    // ������� ������ ����� �����������
    m_additional.clear();

    // ������������� �������������� �����
    setAdditionalDots();

    changeGrid(m_base);
    std::reverse(m_base.begin(), m_base.end());
    changeGrid(m_additional);
    std::reverse(m_additional.begin(), m_additional.end());
}

/**
* �������� ������ ���������� ����������� �� ���������� ������� �����.
* ����� ��������������� ������� ��� ������������ �������� ����� ��� ������������ �����������.
* ������������ ����� ��������� �����������, ����� shiftLinTrigGrid()
*/
BmpVector Grid::intervalMaximums() {
    BmpVector maximums;
    for (auto i = 0; i < m_N_base; ++i) {
        auto last = i == m_N_base - 1 ? m_additional.end() : m_additional.begin() + (i + 1)*m_additionalNCount;
        BmpReal maxElement = *std::max_element(m_additional.begin() + i*m_additionalNCount, last);
        maximums.push_back(maxElement);
    }
    return maximums;
}

