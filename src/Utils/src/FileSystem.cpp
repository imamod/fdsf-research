#include "FileSys.h"

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

    // Считать значения из файла
    json readFile(const std::string& filename) {
        std::ifstream file(filename);
        json jsonData;
        file >> jsonData;
        file.close();
        return jsonData;
    }

    // Вывод значений в файл
    void writeFile(const std::string& filename, const json& data) {
        std::ofstream file(filename);
        file.precision(std::numeric_limits<double>::max_digits10);
        file << data << std::endl;
        file.close();
    }

}
