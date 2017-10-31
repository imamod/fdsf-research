#pragma once

#include "BasicTypes.h"

namespace filesys {
    // �������� ������ � ����
    void writeFile(const std::string& filename, const BmpVector& data);
}

// �������������� ����������
namespace epc {
    typedef BmpReal(*function)(BmpReal x);

    BmpReal Richardson(function f, BmpReal a, BmpReal b, bool countEvery = false);
} // exponential_convergense
