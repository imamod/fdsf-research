#include "BasicService.h"

#include <iostream>
#include <fstream>

namespace filesys {
    // Вывод значений в файл
    void writeFile(const std::string& filename, const BmpVector& data) {
        std::ofstream file(filename);
        file.precision(std::numeric_limits<BmpReal>::max_digits10);
        for (auto const& it : data) {
            file << std::fixed << it << std::endl;
        }
        file.close();
    }
}

namespace test {
    void printVector(const BmpVector& v) {
        for (auto const& it : v) {
            std::cout << it << " ";
        }
        std::cout << std::endl;
    }
}