#include "BasicService.h"

#include <fstream>

namespace filesys {
    // Вывод значений в файл
    void writeFile(const std::string& filename, const BmpVector& data) {
        std::ofstream file(filename);
        file.precision(std::numeric_limits<bmp_real>::max_digits10);
        for (auto const& it : data) {
            file << std::fixed << it << std::endl;
        }
        file.close();
    }
}
