#include "Common.h"
#include "FdHalfQuadratures.h"
#include "Fdsf.h"
#include "FileSys.h"
#include "JsonFields.h"
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

    const double epsilon = 1e-11;
/*
    // TODO: possibly remove
    nlohmann::json calculateWithRichardson(double k, double x) {
        const int N_init = 12;
        int N = N_init;
        double stop_criteria;
        double I_n = quad::euler_maclaurin(x, k, N);
        //print(x, I_n, N);
        // ДЛя проверки по уточнению Ричардсоном
        double diff_prev_InmI2n = 0;
        do {
            double I_2n = quad::euler_maclaurin(x, k, 2 * N);
            stop_criteria = (I_n / I_2n - 1);
            // Сохраняем для уточнению Ричардсоном
            if (N > N_init) {
                double diff = I_2n - I_n;
                double comb = diff / pow(diff_prev_InmI2n, 2);
                std::cout << " N = " << 2 * N << ", diff = " << diff <<
                    ", comb = " << comb << ", check = " << diff*comb / I_2n <<
                    ", check_2 = " << diff*comb << std::endl;
                I_2n *= 1 + comb;
            }
            diff_prev_InmI2n = I_2n - I_n;
            I_n = I_2n;
            N = 2 * N;
            // print(x, I_2n, N);
        } while (abs(stop_criteria) > epsilon);
        nlohmann::json object = nlohmann::json::object();
        object[fd::X] = x;
        object[fd::I] = 2 * I_n;// Смотри формулу (37) препринт 2
        object[fd::N_MAX] = N / 2;
        //std::cout << object.dump() << std::endl;
        return object;
    }
    */

}

TEST_CASE("calculate") {
    nlohmann::json result = nlohmann::json::object();
   // TODO: переработать функцию filesys::createDirectory("quadratures");
    // TODO: setPreciseOutput();
    SECTION("m3half") {
        double x_star = 52;
        calculate(result, fdsf::index::M3_HALF, x_star);
        filesys::writeFile("values_m32.json", result);
    }
    SECTION("mhalf") {
        double x_star = 39;
        calculate(result, fdsf::index::M1_HALF, x_star);
        filesys::writeFile("values_m12.json", result);
    }
    SECTION("half") {
        double x_star = 35;
        calculate(result, fdsf::index::P1_HALF, x_star);
        filesys::writeFile("values_12.json", result);
    }
    SECTION("3half") {
        double x_star = 33;
        calculate(result, fdsf::index::P3_HALF, x_star);
        filesys::writeFile("values_32.json", result);
    }
    SECTION("5half") {
        double x_star = 30;
        calculate(result, fdsf::index::P5_HALF, x_star);
        filesys::writeFile("values_52.json", result);
    }
    SECTION("7half") {
        double x_star = 29;
        calculate(result, fdsf::index::P7_HALF, x_star);
        filesys::writeFile("values_72.json", result);
    }
}

TEST_CASE("richardson_check") {
    nlohmann::json result = nlohmann::json::object();
    SECTION("m3half") {
        double k = -1.5;
        double x_star = 52;
      //  result = calculateWithRichardson(k, x_star);
    }
}