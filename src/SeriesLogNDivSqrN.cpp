#include "SeriesLogNDivSqrN.h"

SeriesLogNDivSqrN::SeriesLogNDivSqrN()
    : m_seriesSum(sum()) {}

SeriesLogNDivSqrN::SeriesLogNDivSqrN(size_t optimal_N)
    : m_N(optimal_N)
    , m_seriesSum(sum()) {}

// ��������� �������� ����� ����
BmpReal SeriesLogNDivSqrN::get() const {
    return m_seriesSum;
}

// ����� ������, �� �������� �������� �����
size_t SeriesLogNDivSqrN::upperBound() const {
    return m_N;
}

// ������ ����� 2<=n<=N-1
BmpReal SeriesLogNDivSqrN::limitedSum() const {
    size_t N = upperBound();
    BmpReal sumValue = log(N) / (2 * pow(N, 2));
    for (size_t n = N - 1; n > 1; --n) {
        sumValue += static_cast<BmpReal>(log(n)) / pow(n, 2);
    }
    return sumValue;
}

/*******************************************************************************
*                                  PRIVATE
*******************************************************************************/

// ������������ �����
BmpReal SeriesLogNDivSqrN::integralPart() const {
    size_t N = upperBound();
    return (1 + log(N)) / N;
};

// ��������� ����� ����
BmpReal SeriesLogNDivSqrN::sum() const {
    return limitedSum() + integralPart() + first() + second();
}

// ������ ������ ��������
BmpReal SeriesLogNDivSqrN::first() const {
    size_t N = upperBound();
    return (2 * log(N) - 1) / (12 * pow(N, 3));
}

// ������ ������ ��������
BmpReal SeriesLogNDivSqrN::second() const {
    size_t N = upperBound();
    return (26 - 24 * log(N)) / (720.0 * pow(N, 5));
}
