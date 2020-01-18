#pragma once

#include "json.hpp"
#include "BasicTypes.h"

using json = nlohmann::json;

namespace filesys {
    // Записать данные в файл
    void writeFile(const std::string& filename, const BmpVector& data);
    // Вывод значений в файл
    void writeFile(const std::string& filename, const std::map<BmpReal, BmpReal>& data);
    // Создать папку
    std::string createDirectory(BmpReal k, size_t n, const std::string& prefix = "");
    // Считать значения из файла
    json readFile(const std::string& filename);
    // Вывод значений в файл
    void writeFile(const std::string& filename, const json& data);
}