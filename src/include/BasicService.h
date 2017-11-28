#pragma once

#include "BasicTypes.h"
#include <functional>

namespace test {
    // ����������� ������
    void printVector(const BmpVector& v, bool printInColumn = false);
    // ����������� �������
    void printMatrix(const BmpMatrix& m);
}

namespace filesys {
    // ������� ������ ������ �� �����
    BmpVector readFile(const std::string& filename);
    // �������� ������ � ����
    void writeFile(const std::string& filename, const BmpVector& data);
    // ������� �����
    std::string createDirectory(BmpReal k, size_t n, const std::string& prefix = "");
}

// �������������� ����������
namespace epc {
    BmpReal Richardson(std::function<BmpReal(const BmpReal&)> f, BmpReal a, BmpReal b, bool countEvery = false);
}

namespace fdsf {
    // ����� �������
    BmpReal NewtonsMethod(BmpReal x, BmpReal k);
}
