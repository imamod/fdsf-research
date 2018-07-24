#pragma once

#include "BasicService.h"
#include <array>

// ��������� �������� ���������
const BmpReal _1 = BmpReal(1);
const BmpReal _2 = BmpReal(2);
const BmpReal _5 = BmpReal(5);
const BmpReal _7 = BmpReal(7);
const BmpReal _691 = BmpReal(691);
const BmpReal _3617 = BmpReal(3617);
const BmpReal _43867 = BmpReal(43867);
const BmpReal _77683 = BmpReal(77683);
const BmpReal _174611 = BmpReal(174611);
const BmpReal _854513 = BmpReal(854513);
const BmpReal _236364091 = BmpReal(236364091);

// �������� pi
inline const BmpReal pi() {
    return boost::math::constants::pi<BmpReal>();
}

// ��������� ������
inline const BmpReal eilerConst() {
    return 0.5772156649015325;
}

// ��������� j
inline const BmpReal j() {
    return 0.76740941382814898;
}

// �������� ����� ���
inline const BmpReal jFromKuzmina() {
    return 0.95024075;
}

// ����� ��������
BmpReal bernulli(int n);

// �����-�������
BmpReal dzetaFunction(int n);
