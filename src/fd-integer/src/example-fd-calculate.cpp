#include "Common.h"
#include "Fdsf.h"
#include "FileSys.h"
#include "JsonFields.h"

#include <functional>

namespace {

    nlohmann::json calculate(std::function<BmpReal(BmpReal)> fd) {
        nlohmann::json params = filesys::readFile("range.json");
        std::cout << params["x_start"] << " " << params["x_end"] << " " << params["span"];
        nlohmann::json result = nlohmann::json::array();
        for (BmpReal x = params["x_start"].get<BmpReal>();
                     x <= params["x_end"].get<BmpReal>();
                     x = x + params["span"].get<BmpReal>()) {
            BmpReal f = fd(x);
            std::cout << "x = " << x << std::endl;
            result.push_back(f);
        }
        return result;
    }
}

TEST_CASE("calc") {
    SECTION("0") {
        INFO("Вычисление диапазона значений функции ФД индекса k = 0");
        auto res = calculate(fdsf::fd_0);
        filesys::writeFile("fd_0.json", res);
    }
    SECTION("1") {
        INFO("Вычисление диапазона значений функции ФД индекса k = 1");
        auto res = calculate(fdsf::fd_1);
        filesys::writeFile("fd_1.json", res);
    }
    SECTION("2") {
        INFO("Вычисление диапазона значений функции ФД индекса k = 2");
        auto res = calculate(fdsf::fd_2);
        filesys::writeFile("fd_2.json", res);
    }
    SECTION("3") {
        INFO("Вычисление диапазона значений функции ФД индекса k = 3");
        auto res = calculate(fdsf::fd_3);
        filesys::writeFile("fd_3.json", res);
    }
    SECTION("4") {
        INFO("Вычисление диапазона значений функции ФД индекса k = 4");
        auto res = calculate(fdsf::fd_4);
        filesys::writeFile("fd_4.json", res);
    }
}