#pragma once

#include "BasicService.h"
#include "Constants.h"

class AsymptoticSeries {
    public:
        AsymptoticSeries(BmpReal k, BmpReal x);

        // �������� ����������� �������� ����� ���� � ����� �
        BmpReal get() const;

    private:
        // ����������� N
        size_t m_N{ 6 };
        // �������� ������������� � ���������������� ����
        const BmpVector m_A = {
            pow(pi(), 2) / 6.0,
            pow(pi(), 4) / 90.0,
            pow(pi(), 6) / 945.0,
            pow(pi(), 8) / 9450.0,
            pow(pi(), 10) / 93555.0,
            691.0 * pow(pi(), 12) / 638512875.0
        };
        // �������� ����
        BmpReal m_seriesSum{ 0 };

        // ���������� ������������ �
        BmpVector getC(BmpReal k, BmpReal x);
        // ��������������
        BmpReal sum(BmpReal k, BmpReal x);

};
