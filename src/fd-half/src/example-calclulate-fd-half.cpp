#include "Common.h"
#include "Fdsf.h"
#include "FileSys.h"
#include "JsonFields.h"

#include <functional>
#include <iostream>

#define RES_FILENAME(name) #name

namespace {

    nlohmann::json calculate(std::function<BmpReal(BmpReal)> fd) {
        nlohmann::json params = filesys::readFile("range.json");
        std::cout << params["x_start"] << " " << params["x_end"] << " " <<params["span"];
        nlohmann::json result = nlohmann::json::array();
        for (BmpReal x = params["x_start"].get<BmpReal>();
                     x <= params["x_end"].get<BmpReal>();
                     x = x + params["span"].get<BmpReal>()) {
            BmpReal f = fd(x);
            std::cout << "x = " << x << std::endl;
            result.push_back(f);
        }
        return result;
        // TODO:
        //filesys::writeFile(RES_FILENAME(fd), result);
    }
}

TEST_CASE("calc") {
    SECTION("m3half") {
        INFO("Вычисление диапазона значений функции ФД индекса k = -3/2");
        nlohmann::json res = calculate(fdsf::fd_m3half);
        filesys::writeFile("fd_m3half.json", res);
    }
    SECTION("m1half") {
        INFO("Вычисление диапазона значений функции ФД индекса k = -1/2");
        auto res = calculate(fdsf::fd_m1half);
        filesys::writeFile("fd_m1half.json", res);
    }
    SECTION("1half") {
        INFO("Вычисление диапазона значений функции ФД индекса k = 1/2");
        auto res = calculate(fdsf::fd_1half);
        filesys::writeFile("fd_1half.json", res);
    }
    SECTION("3half") {
        INFO("Вычисление диапазона значений функции ФД индекса k = 3/2");
        auto res = calculate(fdsf::fd_3half);
        filesys::writeFile("fd_3half.json", res);
    }
    SECTION("5half") {
        INFO("Вычисление диапазона значений функции ФД индекса k = 5/2");
        auto res = calculate(fdsf::fd_5half);
        filesys::writeFile("fd_5half.json", res);
    }
    SECTION("7half") {
        INFO("Вычисление диапазона значений функции ФД индекса k = 7/2");
        auto res = calculate(fdsf::fd_7half);
        filesys::writeFile("fd_7half.json", res);
    }
}
