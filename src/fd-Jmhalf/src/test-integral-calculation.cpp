#include "Common.h"
#include "json.hpp"
#include "JsonFields.h"
#include <chrono>
#include <limits>
#include "Logger.h"
#include "FileSys.h"
#include "EXP_THETA.h"

namespace {
    void print(double x, double I, int N, double stop_criteria = 0.0) {
        std::cout.precision(std::numeric_limits<double>::max_digits10);
        std::cout << "x = " << x << " N = " << N <<
                     " I = " << I << " stop = " << stop_criteria << std::endl;
    }

    std::map<BmpReal, BmpReal> EXP_TAU2_THETA2 = {
    };

    using FermiFunction = std::function<double(double x, double tau, double theta)>;

    // Веса функции
    double weights(int n, int m) {
        //Logger log("weights");
        //log.info("n = " + std::to_string(n) + " m = " + std::to_string(m));
        // в углу вес g = 1/8
        if ((n == m) && !n) {
            return 1.0 / 8;
        }
        // Половинки на диагонали и вертикали g = 1/2
        if ((n == m) || (!n && m)) {
            return 1.0 / 2;
        }
        // Для остальных вес g = 1
        return 1;
    }

    // Подынтегральная функция интегральной ФД
    double fd_Jmhalf(double exp_x, double tau, double theta) {
       // Logger logger("fd_Jmhalf");
        //double a = exp(tau*tau);
        double a = EXP2_TAU_THETA[tau];
        double t = exp_x;
        if (tau == theta) {
            return log(1 + t / a) - t /(a + t) ;
        }
        //double b = exp(theta*theta);
        double b = EXP2_TAU_THETA[theta];
        double denom = a - b;
        double num = a*log(1 + t / a) - b*log(1 + t / b);
        //logger.info("tau = " + std::to_string(tau) + " x = " + std::to_string(x));
        return num / denom;
    }

    // Euler-Macloren Formulas
    double trapz(FermiFunction f, double x, int N) {
        // Верхняя граница, пока непонятно, какая будет
        double T = 12;
        double h = T / N;
        double I = 0.0;
        std::ofstream fout("Jmhalf.txt");
        for (int n = N; n > -1; --n) {
            for (int m = N; m >= n; --m) {
                double tau = n * h;
                double theta = m * h;
             //   EXP_TAU2_THETA2[tau] = exp(tau*tau);
             //   EXP_TAU2_THETA2[theta] = exp(theta*theta);
                I += weights(n, m) * f(x, n * h, m * h);
            }
        }
        // Умножаем на h^2, т.к. двойно интеграл
        // Умножаем на 2 , т.к. интегрируем по треугольнику
        return 2*h*h*I;
    }

    // Критерий останова для формул Эйлера-Маклорена
    const double epsilon = 1e-11;

    nlohmann::json calculate(double x) {
        const int N_init = 12;
        int N = N_init;
        double stop_criteria;
        double exp_x = exp(x);
        double I_n = trapz(fd_Jmhalf, exp_x, N);
        print(x, I_n, N);
        do {
            double I_2n = trapz(fd_Jmhalf, exp_x, 2 * N);
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

    void checkAccuracy(BmpReal x, BmpReal left, BmpReal right) {
        setPreciseOutput();
        std::cout << " x = " << x << ", d = " << left / right - 1 << std::endl;
    }
}

TEST_CASE("calculate") {
    setPreciseOutput();
    std::ofstream f("fd_Jmhalf");
    std::ofstream f_exp("JmExpValues");
    f_exp.precision(std::numeric_limits<BmpReal>::max_digits10);
    nlohmann::json values = nlohmann::json::array();
    size_t i = 0;
    for (double x = 0.0; x <= 60; x+=5) {
    //for (double x = -2.0; x < 1; x++) {
    //    auto start = std::chrono::steady_clock::now();
        values.push_back(calculate(x));
      /*  auto end = std::chrono::steady_clock::now();
        auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "calc x = " << x << " time: " << elapsed_ms.count() << "ms" << std::endl;
        */f << values[i++][fd::I] << std::endl;
    }
    filesys::writeFile("valuesJm2.json", values);
    for (auto x : EXP_TAU2_THETA2) {
        f_exp << "{ " << x.first << ", "<< x.second << " },"<< std::endl;
    }
    f_exp.close();
    f.close();
    //result[fd::RESULT] = values;
}

// Проверка точности вычисления рассчитанных и закешированнных экспонент
TEST_CASE("checkAccuracyExp") {
    const BmpVector cached = {
        0.783238669833194,
        46.58918751923861,
        194.03517728526802,
        442.6570411953065,
        791.6965457165452,
        1240.9561165133946,
        1790.3529161501417,
        2439.8437531948134,
        3189.403142203764,
        4039.014755720472,
        4988.667493982119,
    };

    const BmpVector TRUE_VALUES = {
        0.783238669833194,
        46.58918751923861,
        194.03517728526802,
        442.6570411953065,
        791.6965457165452,
        1240.9561165133946,
        1790.3529161501417,
        2439.8437531948134,
        3189.403142203764,
        4039.014755720472,
        4988.667493982119,
    };
    for (auto lhs = TRUE_VALUES.begin(), rhs = cached.begin();
            lhs != TRUE_VALUES.end() && rhs != cached.end();
            ++lhs, ++rhs) {
        setPreciseOutput();
        std::cout << " d = " << *lhs / *rhs - 1 << std::endl;
    }
}

TEST_CASE("accuracy") {
    const BmpVector GSL_RESULT = {
        0.78323866983319235,
        46.589187519238614,
        194.03517728526793,
        442.65704119530864,
        791.69654571653882,
        1240.9561165133962,
        1790.352916150179,
        2439.843753194812,
        3189.4031422038,
        4039.014755720504,
        4988.667493982043,
    };
    double x = 0.0;
    for (auto gslValue : GSL_RESULT) {
        nlohmann::json result = calculate(x);
        checkAccuracy(x, gslValue, result[fd::I]);
        x += 5;
    }
}

namespace {
    // Кэшированные значения экспонент для tau и theta.
    std::map<double, double> exp2Values;
    // Кэшированные значения логарифмов log(1 + t / x), где x = {a, b}
    std::map<double,double> logValues;

    /**
     * Оптимизированный метод вычисления J(x).
     * in: exp(x)
     * in: точка по tau
     * in: точка по theta
     */
    double fd_Jmhalf_modified(double exp_x, double tau, double theta) {
        // Если нет такого значения в кэшированном, вычисляем и сохраняем
        if (!exp2Values[tau]) {
            exp2Values[tau] = exp(tau*tau);
        }
        double a = exp2Values[tau];
        double t = exp_x;
        if (tau == theta) {
            // Если log для данного a не рассчитан, вычисляем и сохраняем
            /*if (!logValues[a]) {
                logValues[a] = log(1 + t / a);
            }
            return logValues[a] - t / (a + t);*/
            return log(1 + t / a) - t / (a + t);
        }
        // Если нет такого значения в кэшированном, вычисляем и сохраняем
        // Расчет ведется на квадратной равномерной сетке, поэтому диапазон значений tau и theta одинаков
        if (!exp2Values[theta]) {
            exp2Values[theta] = exp(theta*theta);
        }
        double b = exp2Values[theta];
        // Если log для данного b не рассчитан, вычисляем и сохраняем
       /* if (!logValues[b]) {
            logValues[b] = log(1 + t / b);
        }
        double num = a*logValues[a] - b*logValues[b];*/
        double num = a*log(1 + t / a) - b*log(1 + t / b);
        return num / (a - b);
    }

    double modifiedTrapz(FermiFunction f, double exp_x, int N) {
        // Верхняя граница, пока непонятно, какая будет
        double T = 12;
        double h = T / N;
        double I = 0.0;
        for (int n = N; n > -1; --n) {
            std::vector<double> u;
            for (int m = N; m >= n; --m) {
                u.push_back(weights(n, m) * f(exp_x, n * h, m * h));
            }
            for (auto it = u.begin(); it != u.end(); ++it) {
                I += *it;
            }
        }
        // Умножаем на h^2, т.к. двойной интеграл
        // Умножаем на 2 , т.к. интегрируем по треугольнику
        return 2 * h*h*I;
    }
}

TEST_CASE("check") {
    setPreciseOutput();
    auto begin = std::chrono::steady_clock::now();
    for (auto x : {0, 5}) {
        double exp_x = exp(x);
        double I = trapz(fd_Jmhalf, exp_x, 12);
        std::cout << I << std::endl;
    }
    auto end = std::chrono::steady_clock::now();
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    std::cout << "non modifyed trapz = " << elapsed_ms.count() << std::endl;

    begin = std::chrono::steady_clock::now();
    for (auto x : { 0, 5 }) {
        double exp_x = exp(x);
        double I = modifiedTrapz(fd_Jmhalf_modified, exp_x, 12);
        std::cout << "x = " << x << ", I = " << 4*I << std::endl;
    }
    end = std::chrono::steady_clock::now();
    elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    std::cout << "modifyed trapz = " << elapsed_ms.count() << std::endl;
}