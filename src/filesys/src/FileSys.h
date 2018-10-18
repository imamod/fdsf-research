#pragma once

#include "BasicTypes.h"

namespace filesys {
    // ������� ������ ������ �� �����
    BmpVector readFile(const std::string& filename);
    // �������� ������ � ����
    void writeFile(const std::string& filename, const BmpVector& data);
    // ������� �����
    std::string createDirectory(BmpReal k, size_t n, const std::string& prefix = "");
}