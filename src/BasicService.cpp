#include "BasicService.h"

#include <iostream>
#include <fstream>
#include <experimental/filesystem>


namespace filesys {
    // Считать вектор данных из файла
    BmpVector readFile(const std::string& filename) {
        std::ifstream file(filename, std::ios::binary);
        file.precision(std::numeric_limits<BmpReal>::max_digits10);
        BmpVector data;
        while (true) {
            if (file.eof()) {
                break;
            }
            BmpReal value;
            file >> value;
            data.push_back(value);
        }
        file.close();
        return data;
    }

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
#ifdef HIGH_PRECISION
        std::string path(prefix + "k");
#else
        std::string path(prefix + "k" + std::to_string(k));
#endif // HIGH_PRECISION
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