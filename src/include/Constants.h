#pragma once

#include "BasicService.h"

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
