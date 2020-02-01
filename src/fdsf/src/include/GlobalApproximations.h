#pragma once

#include "BasicTypes.h"

enum GlobalFive {
    ONE_COEFFICIENT = 1,
    LOW_TEMP,
    BEST_PREC
};

namespace global_formula {
    /**
     * ��������� ������� �� �� ���������� 5-������� ��������
     * @param ����� x, � ������� ����� ��������� �������� ������� ��
     * @param ������ ������� ��
     * @param ��� ���������� �������
     * @return ����������� �������� ������� ��
     */
    BmpReal calculateFive(BmpReal x, BmpReal k, uint8_t approxType);
}
