#include "BasicService.h"

#include <iostream>
#include <fstream>
#include <experimental/filesystem>


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

    std::string createDirectory(BmpReal k, size_t n, const std::string& prefix) {
        namespace fs = std::experimental::filesystem;
        if (!prefix.empty() && !fs::exists(prefix)) {
            fs::create_directory(prefix);
        }
        std::string path(prefix + "k" + std::to_string(k));
        if (!fs::exists(path)) {
            fs::create_directory(path);
        }
        path += "/n" + std::to_string(n);
        if (!fs::exists(path)) {
            fs::create_directory(path);
        }
        return path + "/";
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