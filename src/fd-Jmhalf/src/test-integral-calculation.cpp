#include "Common.h"
#include "json.hpp"
#include "JsonFields.h"
#include "Logger.h"

namespace {
    void print(double x, double I, int N, double stop_criteria = 0.0) {
        std::cout.precision(std::numeric_limits<double>::max_digits10);
        std::cout << "x = " << x << " N = " << N <<
                     " I = " << I << " stop = " << stop_criteria << std::endl;
    }

    using FermiFunction = std::function<double(double x, double tau, double theta)>;

    // Веса функции
    double weights(int n, int m) {
        Logger log("weights");
        log.info("n = " + std::to_string(n) + " m = " + std::to_string(m));
        // в углу вес g = 1/8
        if ((n == m) && !n) {
            return 1.0 / 8;
        }
        // Половинки на диагонали и вертикали g = 1/2
        if ((n == m) ||
            (!n && m)) {
            return 1.0 / 2;
        }
        // Для остальных вес g = 1
        return 1;
    }

    // Подынтегральная функция интегральной ФД
    double fd_Jmhalf(double x, double tau, double theta) {
        Logger logger("fd_Jmhalf");
        double a = exp(tau*tau);
        double t = exp(x);
        if (tau == theta) {
            return log(1 + t / a) - t /(a + t) ;
        }
        double b = exp(theta*theta);
        double denom = a - b;
        double num = a*log(1 + t / a) - b*log(1 + t / b);
        //log.info("tau = " + std::to_string(tau) + " x = " + std::to_string(x) + " ch(x) = " + std::to_string(ch_x));
        return num / denom;
    }

    // Euler-Macloren Formulas
    double trapz(FermiFunction f, double x, double T, size_t N) {
        Logger log("trapz");
        double h = T / N;
        double I = 0.0;
        for (int n = N; n > -1; --n) {
            std::vector<double> u;
            for (int m = N; m >= n; --m) {
                u.push_back(weights(n, m) * f(x, n * h, m * h));
                log.info(std::to_string(u.back()));
            }
            // weights(n , m) *
            for (auto it = u.begin(); it != u.end(); ++it) {
                I += *it;
            }
        }
        // Умножаем на h^2, т.к. двойно интеграл
        // Умножаем на 2 , т.к. интегрируем по треугольнику
        return 2*h*h*I;
    }

    double euler_maclaurin(double x, int N, double T) {
        FermiFunction f = fd_Jmhalf;
        return trapz(f, x, T, N);
    }

    // Критерий останова для формул Эйлера-Маклорена
    const double epsilon = 1e-11;

    nlohmann::json calculate(double x) {
        const int N_init = 12;
        int N = N_init;
        double stop_criteria;
        // Верхняя граница, пока непонятно, какая будет
        // TODO: Удалить
        double T = 12;
        double I_n = euler_maclaurin(x, N, T);
        print(x, I_n, N);
        do {
            double I_2n = euler_maclaurin(x, 2 * N, T);
            stop_criteria = (I_n / I_2n - 1);
            I_n = I_2n;
            N = 2 * N;
            print(x, I_2n, N, stop_criteria);
        } while (abs(stop_criteria) > epsilon);
        nlohmann::json object = nlohmann::json::object();
        object[fd::X] = x;
        object[fd::I] = 4 * I_n;// Смотри формулу (37) препринт 2
        object[fd::N_MAX] = N / 2;
        //std::cout << object.dump() << std::endl;
        return object;
    }
}

TEST_CASE("calculate") {
    std::ofstream f("fd_Jmhalf");
    nlohmann::json values = nlohmann::json::array();
    for (double x = 0.0; x <= 50; x+=5) {
        values.push_back(calculate(x));
        f << values[x][fd::I] << std::endl;
    }
    f.close();
    //result[fd::RESULT] = values;
}