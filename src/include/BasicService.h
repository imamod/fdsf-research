#pragma once

#include "BasicTypes.h"

namespace filesys {
    // �������� ������ � ����
    void writeFile(const std::string& filename, const BmpVector& data);
}

// �������������� ����������
namespace epc {
    typedef bmp_real(*function)(bmp_real x);

    bmp_real Richardson(function f, bmp_real a, bmp_real b, bool countEvery = false);
} // exponential_convergense
