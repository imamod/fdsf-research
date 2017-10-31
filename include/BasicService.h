#pragma once

#include "BasicTypes.h"

namespace filesys {
    // Записать данные в файл
    void writeFile(const std::string& filename, const BmpVector& data);
}
