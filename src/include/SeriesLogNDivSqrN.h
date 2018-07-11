#pragma once

#include "BasicTypes.h"

// �ODO: class series ���
/**
* ����� ������� ���� log(n)/(n^2), ��������� ��� ���������� �������������� ��� x>>1
*/
class SeriesLogNDivSqrN {
    public:
        SeriesLogNDivSqrN();
        SeriesLogNDivSqrN(size_t optimal_N);

        // ��������� ��������
        BmpReal get() const;
        // ����� ������, �� �������� �������� �����
        size_t upperBound() const;
        // ������� ����� 2<=n<=N-1
        BmpReal limitedSum() const;

    private:
        // ����������� N
        size_t m_N{ 300 };
        // �������� ����
        BmpReal m_seriesSum{ 0 };

        // ���������� ����� ����
        BmpReal sum() const;
        // ������ ������ ��������
        BmpReal first() const;
        // ������ ������ ��������
        BmpReal second() const;
        // ������������ �����
        BmpReal integralPart() const;
};
