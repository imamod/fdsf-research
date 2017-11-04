#pragma once

#include "BasicTypes.h"

class Grid {
    public:
        Grid(size_t N, size_t addNCount = 11);

        // ���������� ����� ������� ����� �� ��������� ������
        void setLinearGrid();

        // ���������� ����� ������� ����� �������-������������������� ������ �����
        void setLinearTrigonometricGrid();

        /**
         * ������ ������ �������-������������������ ����� � ������� �����, ���� 10(?)
         * �������������� ����� ����� ������ ����� ������� �����. ��������� ������ ��� ��������� ��������
         */
        void setLinearTrigonometricGridRight();

        // �������� ������ ������� �����
        BmpVector base() const;

        // �������� ������ �������������� ����� (������ � ��������)
        BmpVector additional() const;

        // �������� x �� y
        BmpVector xByY(const BmpVector& y);

    private:
        // ����� ������� �����
        size_t m_N_base;
        // ����� �������������� ����� ����� ��������
        size_t m_additionalNCount;
        // ������ ������� �����
        BmpVector m_base;
        // ������ �������������� �����
        BmpVector m_additional;
};
