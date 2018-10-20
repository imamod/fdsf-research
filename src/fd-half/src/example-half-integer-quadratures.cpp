#include "Common.h"
#include "json.hpp"
#include "FileSys.h"

#include <fstream>
#include <iostream>
#include <functional>

using json = nlohmann::json;

namespace fd {
    const std::string X = "x";
    const std::string I = "I";
    const std::string K = "k";
    const std::string N_MAX = "N";
    const std::string RESULT = "result";
}

namespace {

    void print(double x, double I, int N) {
        std::cout.precision(std::numeric_limits<double>::max_digits10);
        std::cout << "x = " << x << " N = " << N  << " I = "<< I << std::endl;
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

    using FermiFunction = std::function<double(const double& x, double k, double tau)>;

    // Подынтегральная функция для индекса k > -3/2
    double fd_m12(double tau, double x, double k) {
        double denom = 1 + exp(tau * tau - x);
        return pow(tau, 2*k + 1) / denom;
    }

    /*
    // Подынтегральная функция для индекса k = -3/2
    double fermi_dirak_m3half( double x,
        double k) {
        double exp_ksi = exp(-a * ksi * ksi / (1 - ksi * ksi));
        double sum_exp = exp_ksi + exp(-x);
        double exp_diff = exp(-a * ksi * ksi / (1 - ksi * ksi) - x);

        return (-4 * exp_diff * sqrt(abs(a)) * exp_ksi) /
            (pow(1 - ksi * ksi, double(3.0 / 2)) * sum_exp * sum_exp);
    }
    */
    // Euler-Macloren Formulas
    double trapz(FermiFunction f, double x, const double k, size_t N) {
        double h = 12.0 / N;
        std::vector<double> u;
        for (size_t i = 0; i < N + 1; ++i) {
            u.push_back(f(i * h, x, k));
        }
        u.at(0) /= 2.0;
        u.at(N) /= 2.0;
        double I = 0.0;
        for (auto it = u.rbegin(); it != u.rend(); ++it) {
            I += *it;
        }
        //double I = std::accumulate(u.rbegin(), u.rend(), 0.0);
        return h*I;
    }

    double euler_maclaurin(double x, double k, int N) {
        /*FermiFunction f = (k == -3.0 / 2) ? fermi_dirak_m3half :
            fermi_dirak_half_integer; */
        return trapz(fd_m12, x, k, N);
    }

    json calculate(double k, double x) {
        int N = 12;
        double stop_criteria;
        double I_n = euler_maclaurin(x, k, N);
        print(x, I_n, N);
        do {
            double I_2n = euler_maclaurin(x, k, 2 * N);
            stop_criteria = (I_n / I_2n - 1);
            I_n = I_2n;
            N = 2 * N;
            print(x, I_2n, N);
        } while (abs(stop_criteria) > 1e-11);
        json object = json::object();
        object[fd::X] = x;
        object[fd::I] = 2*I_n;// Смотри формулу (37) препринт 2
        object[fd::N_MAX] = N/2;
        //std::cout << object.dump() << std::endl;
        return object;
    }

    void calculate(json& result, double k, double x_star) {
        std::ofstream f(std::to_string(k));
        json values = json::array();
        for (double x = 0.0; x < x_star + 1; ++x) {
            values.push_back(calculate(k, x));
            f << values[x][fd::I] << std::endl;
        }
        f.close();
        result[fd::K] = k;
        result[fd::RESULT] = values;
    }
}

TEST_CASE("half") {
    json result = json::object();
   // TODO: переработать функцию filesys::createDirectory("quadratures");
    // TODO: setPreciseOutput();
    SECTION("m12") {
        double k = -0.5;
        double x_star = 39;
        calculate(result, k, x_star);
        writeFile("values_m12.json", result);
    }
    SECTION("12") {
        double k = 0.5;
        double x_star = 35;
        calculate(result, k, x_star);
        writeFile("values_12.json", result);
    }
    SECTION("32") {
        double k = 1.5;
        double x_star = 33;
        calculate(result, k, x_star);
        writeFile("values_32.json", result);
    }
    SECTION("52") {
        double k = 2.5;
        double x_star = 30;
        calculate(result, k, x_star);
        writeFile("values_52.json", result);
    }
    SECTION("72") {
        double k = 3.5;
        double x_star = 29;
        calculate(result, k, x_star);
        writeFile("values_72.json", result);
    }
}