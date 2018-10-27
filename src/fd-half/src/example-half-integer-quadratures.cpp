#include "Common.h"
#include "FdHalfQuadratures.h"
#include "FileSys.h"
#include <fstream>

namespace {

    void calculate(nlohmann::json& result, double k, double x_star) {
        std::ofstream f(std::to_string(k));
        nlohmann::json values = nlohmann::json::array();
        for (double x = 0.0; x < x_star + 1; ++x) {
            values.push_back(quad::calculate(k, x));
            f << values[x][fd::I] << std::endl;
        }
        f.close();
        result[fd::K] = k;
        result[fd::RESULT] = values;
    }
}

TEST_CASE("calculate") {
    nlohmann::json result = nlohmann::json::object();
   // TODO: переработать функцию filesys::createDirectory("quadratures");
    // TODO: setPreciseOutput();
    SECTION("m3half") {
        double k = -1.5;
        double x_star = 52;
        calculate(result, k, x_star);
        filesys::writeFile("values_m32.json", result);
    }
    SECTION("mhalf") {
        double k = -0.5;
        double x_star = 39;
        calculate(result, k, x_star);
        filesys::writeFile("values_m12.json", result);
    }
    SECTION("half") {
        double k = 0.5;
        double x_star = 35;
        calculate(result, k, x_star);
        filesys::writeFile("values_12.json", result);
    }
    SECTION("3half") {
        double k = 1.5;
        double x_star = 33;
        calculate(result, k, x_star);
        filesys::writeFile("values_32.json", result);
    }
    SECTION("5half") {
        double k = 2.5;
        double x_star = 30;
        calculate(result, k, x_star);
        filesys::writeFile("values_52.json", result);
    }
    SECTION("7half") {
        double k = 3.5;
        double x_star = 29;
        calculate(result, k, x_star);
        filesys::writeFile("values_72.json", result);
    }
}