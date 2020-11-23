#include "Common.h"
#include "FdHalfQuadratures.h"
#include "FdIndex.h"
#include "FileSys.h"
#include "JsonFields.h"
#include <fstream>

#include "Richardson.h"
#include "FermiDirakFunction.h"

namespace {

    // Вычисление констант exp(tau*tau)
    void calculateTau2() {
        std::map<BmpReal, BmpReal> result;
        size_t N = 12;
        double h = 12.0 / N;
        for (int i = N - 1; i > 0; --i) {
            double currentTau = i*h;
            result[currentTau] = exp(currentTau*currentTau);
        }
        do {
            N *= 2;
            h = 12.0 / N;
            for (int i = N - 1; i > 0; i=i-2) {
                double currentTau = i*h;
                result[currentTau] = exp(currentTau*currentTau);
            }
        } while ( N != 2*768);
        filesys::writeFile("constTau.json", result);
    }


    void calculate(nlohmann::json& result, double k, double x_star) {
        std::ofstream f(std::to_string(k));
        nlohmann::json values = nlohmann::json::array();
        // Расчет реперной точки x = -2
        values.push_back(quad::calculate(k, -2));
        f << values[0][fd::I] << std::endl;
        // Расчет реперной точки x = -1
        values.push_back(quad::calculate(k, -1));
        f << values[1][fd::I] << std::endl;
        for (double x = 0.0; x < x_star + 1; x = x + 5) {
            values.push_back(quad::calculate(k, x));
            f << values[x + 2][fd::I] << std::endl;
        }
        f.close();
        result[fd::K] = k;
        result[fd::RESULT] = values;
    }

    const double epsilon = 1e-11;

    bool isOldCriteria(double prev, double cur) {
        return abs(prev/cur - 1) < epsilon;
    }

    bool isNewCriteria(double pprev, double prev, double cur) {
        double lhs = cur - prev;
        double rhs = prev - pprev;
        double res = pow(lhs, 3) / (cur * pow(rhs, 2));
        return (res < epsilon) && (lhs*rhs) > 0;
    }

    // TODO: possibly remove
    nlohmann::json calculateWithRichardson(double k, double x) {
        const int N_init = 12;
        const double epsilon = 1e-11;
        RichardsonResult m_result{ N_init, 0 };
        FermiDirakFunction m_func(x, k);
        m_result.I = m_func.calculate(m_result.N);
        std::cout << "N = " << m_result.N << ", I = " << m_result.I << std::endl;
        double I_pprev;
        for (int count = 0;; ++count) {
            m_result.N *= 2;
            double I_2n = m_func.calculate(m_result.N, m_result.I);
            if (isOldCriteria(m_result.I, I_2n)) {
                m_result.I = I_2n;
                std::cout << "Old"<< std::endl;
                break;
            }
            if (count > 1 && isNewCriteria(I_pprev, m_result.I, I_2n)) {
                m_result.I = I_2n + (pow((I_2n - m_result.I), 3)) / (pow((m_result.I - I_pprev), 2));
                std::cout << "New" << std::endl;
                break;
            }
            I_pprev = m_result.I;
            m_result.I = I_2n;
            std::cout << "N = " << m_result.N << ", I = " << m_result.I << std::endl;
        }
        // Домножаем значение интеграла на коэффициент перед ним ( смотри формулы (30, 34) препринт 2 )
        const BmpReal coeff = (fdsf::index::M3_HALF == m_func.index()) ? -1 : 2;
        m_result.I *= coeff;
        std::cout << "final I = " << m_result.I << std::endl;
        nlohmann::json object = nlohmann::json::object();
        object[fd::X] = x;
        object[fd::I] = m_result.I;
        return object;
    }

}

TEST_CASE("calculate") {
    nlohmann::json result = nlohmann::json::object();
   // TODO: переработать функцию filesys::createDirectory("quadratures");
    // TODO: setPreciseOutput();
    double x_star = 60;
    SECTION("m3half") {
        //double x_star = 52;
        calculate(result, fdsf::index::M3_HALF, x_star);
        filesys::writeFile("values_m32.json", result);
    }
    SECTION("mhalf") {
        //double x_star = 39;
        calculate(result, fdsf::index::M1_HALF, x_star);
        filesys::writeFile("values_m12.json", result);
    }
    SECTION("half") {
        //double x_star = 35;
        calculate(result, fdsf::index::P1_HALF, x_star);
        filesys::writeFile("values_12.json", result);
    }
    SECTION("3half") {
        //double x_star = 33;
        calculate(result, fdsf::index::P3_HALF, x_star);
        filesys::writeFile("values_32.json", result);
    }
    SECTION("5half") {
        //double x_star = 30;
        calculate(result, fdsf::index::P5_HALF, x_star);
        filesys::writeFile("values_52.json", result);
    }
    SECTION("7half") {
        //double x_star = 29;
        calculate(result, fdsf::index::P7_HALF, x_star);
        filesys::writeFile("values_72.json", result);
    }
}

TEST_CASE("new_stop_criteria") {
    setPreciseOutput();
    auto OLD_CRITERIA = filesys::readFile("values_12.json");
    auto packs = OLD_CRITERIA[fd::RESULT];
    for (auto x : {0.0, 5.0, 10.0, 15.0, 20.0, 25.0}) {
        std::cout << "x = " << x << std::endl;
        nlohmann::json result = calculateWithRichardson(fdsf::index::P1_HALF, x);
        std::cout << result.dump() << std::endl;
        for (auto pack : packs) {
            if (pack[fd::X] == x) {
                double oldValue = pack[fd::I];
                double newValue = result[fd::I];
                std::cout << "relative error = " << abs(newValue / oldValue - 1) << std::endl;
                break;
            }
        }
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

TEST_CASE("constTau") {
    calculateTau2();
}
