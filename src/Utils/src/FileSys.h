#pragma once

#include "json.hpp"
#include "BasicTypes.h"

using json = nlohmann::json;

namespace filesys {
    // �������� ������ � ����
    void writeFile(const std::string& filename, const BmpVector& data);
    // ������� �����
    std::string createDirectory(BmpReal k, size_t n, const std::string& prefix = "");
    // ������� �������� �� �����
    json readFile(const std::string& filename);
    // ����� �������� � ����
    void writeFile(const std::string& filename, const json& data);
}