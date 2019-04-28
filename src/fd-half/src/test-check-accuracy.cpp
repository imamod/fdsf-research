#include "Common.h"
#include "FileSys.h"

namespace {

    const BmpReal X_LEFT = 0;

    void checkAccuracy(BmpReal k, BmpReal x, BmpReal left, BmpReal right) {
        setPreciseOutput();
        std::cout << "k = " << k << ", x = " << x << ", d = " << left/right - 1 << std::endl;
    }
}

TEST_CASE("check_fcs_and_quad") {
    SECTION("m_3half") {
        BmpReal I_fcs = -1.3474364777155079;
        BmpReal I_quad = -1.3474364777155081;
        BmpReal k = -1.5;
        checkAccuracy(k, X_LEFT, I_fcs, I_quad);
    }
    SECTION("m_half") {
        BmpReal I_fcs = 1.0721549299401913;
        BmpReal I_quad = 1.0721549299401913;
        BmpReal k = -0.5;
        checkAccuracy(k, X_LEFT, I_fcs, I_quad);
    }
    SECTION("half") {
        BmpReal I_fcs = 0.67809389515310081;
        BmpReal I_quad = 0.6780938951531011;
        BmpReal k = 0.5;
        checkAccuracy(k, X_LEFT, I_fcs, I_quad);
    }
    SECTION("3half") {
        BmpReal I_fcs =  1.1528038370883613;
        BmpReal I_quad = 1.1528038370883615;
        BmpReal k = 1.5;
        checkAccuracy(k, X_LEFT, I_fcs, I_quad);
    }
    SECTION("5half") {
        BmpReal I_fcs =  3.0825860828374179;
        BmpReal I_quad = 3.0825860828374188;
        BmpReal k = 2.5;
        checkAccuracy(k, X_LEFT, I_fcs, I_quad);
    }
    SECTION("7half") {
        BmpReal I_fcs =  11.183716751693318;
        BmpReal I_quad = 11.18371675169332;
        BmpReal k = 3.5;
        checkAccuracy(k, X_LEFT, I_fcs, I_quad);
    }
}

TEST_CASE("check_asym_and_quad") {
    SECTION("m_3half") {
        // x = 44
        BmpReal I_quad = -0.30170449383410713;
        //BmpReal I_asym = -0.301704493834138; //N =10
        BmpReal I_asym = -0.30170449383410813; // N = 11
        checkAccuracy(-1.5, 44, I_asym, I_quad);
        // x = 52
        I_quad = -0.27747711522431084;
        I_asym = -0.27747711522431112;
        checkAccuracy(-1.5, 52, I_asym, I_quad);
    }
    SECTION("m_half") {
        BmpReal I_quad = 12.486609377850819;
        BmpReal I_asym = 12.486609377850817; // N = 9
        //BmpReal I_asym = 12.486609377850813; // N = 10
        checkAccuracy(-0.5, 39, I_asym, I_quad);
    }
    SECTION("half") {
        BmpReal I_quad = 138.18098265902458;
        BmpReal I_asym = 138.18098265902452;
        BmpReal k = 0.5;
        checkAccuracy(k, 35, I_asym, I_quad);
    }
    SECTION("3half") {
        BmpReal I_quad = 2516.501868640833;
        BmpReal I_asym = 2516.501868640833; // N = 9
        //BmpReal I_asym = 2516.5018686408321; // N = 10
        BmpReal k = 1.5;
        checkAccuracy(k, 33, I_asym, I_quad);
    }
    SECTION("5half") {
        BmpReal I_quad = 42929.25758509996;
        BmpReal I_asym = 42929.257585099935;
        BmpReal k = 2.5;
        checkAccuracy(k, 30, I_asym, I_quad);
    }
    SECTION("7half") {
        BmpReal I_quad = 872613.5640130762;
        BmpReal I_asym = 872613.56401307601;
        BmpReal k = 3.5;
        checkAccuracy(k, 29, I_asym, I_quad);
    }
}

/**
 * Обоснование замены вычисления экспоненты вычислением корня 4 степени
 */
namespace {

    /* Вычиление прямых экспонент */
    std::vector<BmpReal> direct_exponents(size_t N, BmpReal x, BmpReal k) {
        BmpReal T_left = 0, T_right = 12.0;
        BmpReal step = T_right / N;
        BmpVector exponents;
        for (BmpReal t = T_left; t <= T_right; t += step) {
            exponents.push_back(exp(t*t));
        }
        return exponents;
    }

    /* Замена вычисления экспонент квадратным корнем */
    std::vector<BmpReal> sqrtRefinement(size_t N_max, BmpReal x, BmpReal k) {
        size_t N = 1;
        BmpReal T_left = 0, T_right = 12.0;
        // Поправка c = exp(-144) ^ (1/4)
        BmpReal correction = exp( -144.0 / 4);
        std::vector<BmpReal> exponents = { exp(T_left*T_left), exp(T_right*T_right) };
        while (N < N_max) {
            std::vector<BmpReal> values;
            for (int n = 0; n < N; ++n) {
                BmpReal newElement = sqrt(exponents[n] * exponents[n + 1]) * correction;
                values.push_back(newElement);
            }
            exponents.insert(exponents.begin() + exponents.size(), values.begin(), values.end());
            std::sort(exponents.begin(), exponents.end());
            N *= 2;
            correction = pow(correction, 1.0 / 4);
        }
        return exponents;
    }
}

TEST_CASE("accuracyExpVsSqrt") {
    size_t N_max = 1024;
    BmpReal k = 1.0 / 2;
    BmpReal x = 0;
    nlohmann::json directExp = direct_exponents(N_max, x, k);
    filesys::writeFile("accuracy_direct_exp.json", directExp);
    nlohmann::json directExpSqrt = sqrtRefinement(N_max, x, k);
    filesys::writeFile("accuracy_sqrt_exp.json", directExpSqrt);
    nlohmann::json delta;
    for (int i = 0; i < N_max; ++i) {
        delta.push_back(directExp[i].get<BmpReal>() / directExpSqrt[i].get<BmpReal>() - 1);
    }
    filesys::writeFile("accuracy_sqrt_exp_delta.json", delta);
}

#include <chrono>
namespace {

    inline long double sum_check(long double x) {
        return x + x;
    }

    inline long double mul_check(long double x) {
        return x * x;
    }

    long double div_check(long double x) {
        return x / x;
    }

    long double log_check(long double x) {
        return log(x);
    }

    long double exp_check(long double x) {
        return exp(x);
    }

    long double sqrt_check(long double x) {
        return sqrt(x);
    }

    long double sin_check(long double x) {
        return sin(x);
    }

    long double asin_check(long double x) {
        return asin(x);
    }

    long double pow_check(long double x) {
        return pow(x, 0.25);
    }

    const size_t OPERATION_TIMES = 100000000;
    const double SMALL_NUMBER = pi() * exp(1);
    const double LARGE_NUMBER = OPERATION_TIMES * pi() * exp(1);

    void check_epty_cycle() {
        auto start = std::chrono::system_clock::now();
        for (size_t i = 0; i < OPERATION_TIMES; ++i) {}
        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> diff = (end - start);
        std::cout << std::fixed << std::setprecision(16) <<
                    " empty cycle time: " << 1000000000 * diff.count() / OPERATION_TIMES << " ns." << std::endl;
    }

    /* Проверка трудоемкости операции возведения в целую степень */
    void check_integer_pow(double number, size_t powNumber) {
        auto start = std::chrono::system_clock::now();
        double res;
        for (size_t i = 0; i < OPERATION_TIMES; ++i) {
            res = pow(number, powNumber);
        }
        std::cout << res << std::endl;
        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> diff = (end - start);
        std::cout << std::fixed << std::setprecision(16) <<
            " pow " << powNumber << " time: " << 1000000000 * diff.count() / OPERATION_TIMES << " ns." << std::endl;
    }

    /* Проверка трудоемкости операций */
    void checkOperation(double number, const std::string& operationName, std::function<long double(long double)> operation) {
        auto start = std::chrono::system_clock::now();
        for (size_t i = 0; i < OPERATION_TIMES; ++i) {
            operation(number);
        }
        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> diff = (end - start);
        std::cout << std::fixed << std::setprecision(16) << operationName
                  << " time: " << 1000000000 * diff.count() / OPERATION_TIMES << " ns." << std::endl;
    }
}

/**
small number sum time: 179.46 ns.
large number sum time: 179.59 ns.
small number mul time: 179.40 ns.
large number mul time: 180.15 ns.
small number div time: 179.24 ns.
large number div time: 180.50 ns.
small number log time: 232.68 ns.
large number log time: 234.74 ns.
small number exp time: 233.36 ns.
large number exp time: 606.00 ns.
small number sqrt time: 228.72 ns.
large number sqrt time: 227.99 ns.
small number sin time: 237.08 ns.
large number sin time: 271.40 ns.
small number arcsin time: 570.73 ns.
large number arcsin time: 569.29 ns.

pow
empty cycle time: 3.83 ns.
pow 2 time: 147.43 ns.
pow 2 time: 145.54 ns.
pow 3 time: 144.99 ns.
pow 3 time: 145.32ns.
pow 4 time: 145.62 ns.
pow 4 time: 145.28 ns.
pow 5 time: 145.81 ns.
pow 5 time: 145.16 ns.
pow 11 time: 145.35 ns.
pow 11 time: 145.26 ns.
pow 8 time: 145.63 ns.
pow 853973422 time: 403.12 ns.
*/

TEST_CASE("checkAdditional") {
    {
        INFO("Проверка пустого цикла");
        check_epty_cycle();
    }
    {
        INFO("Проверка возведения во 2 степень");
        check_integer_pow(SMALL_NUMBER, 2);
        check_integer_pow(LARGE_NUMBER, 2);
    }
    {
        INFO("Проверка возведения в 3 степень");
        check_integer_pow(SMALL_NUMBER, 3);
        check_integer_pow(LARGE_NUMBER, 3);
    }
    {
        INFO("Проверка возведения в 4 степень");
        check_integer_pow(SMALL_NUMBER, 4);
        check_integer_pow(LARGE_NUMBER, 4);
    }
    {
        INFO("Проверка возведения в 5 степень");
        check_integer_pow(SMALL_NUMBER, 5);
        check_integer_pow(LARGE_NUMBER, 5);
    }
    {
        INFO("Проверка возведения в 11 степень");
        check_integer_pow(SMALL_NUMBER, 11);
        check_integer_pow(LARGE_NUMBER, 11);
    }
    {
        INFO("Проверка возведения в степень");
        check_integer_pow(SMALL_NUMBER, SMALL_NUMBER);
        check_integer_pow(LARGE_NUMBER, LARGE_NUMBER);
    }
}

TEST_CASE("checkOperationTime") {
    {
        INFO("Проверка единичного сложения");
        checkOperation(SMALL_NUMBER, "small number sum", sum_check);
        checkOperation(LARGE_NUMBER, "large number sum", sum_check);
    }
    {
        INFO("Проверка единичного умножения");
        checkOperation(SMALL_NUMBER, "small number mul", mul_check);
        checkOperation(LARGE_NUMBER, "large number mul", mul_check);
    }
    {
        INFO("Проверка единичного деления");
        checkOperation(SMALL_NUMBER, "small number div", div_check);
        checkOperation(LARGE_NUMBER, "large number div", div_check);
    }
    {
        INFO("Проверка log");
        checkOperation(SMALL_NUMBER, "small number log", log_check);
        checkOperation(LARGE_NUMBER, "large number log", log_check);
    }
    {
        INFO("Проверка exp");
        checkOperation(SMALL_NUMBER, "small number exp", exp_check);
        checkOperation(LARGE_NUMBER, "large number exp", exp_check);
    }
    {
        INFO("Проверка sqrt");
        checkOperation(SMALL_NUMBER, "small number sqrt", sqrt_check);
        checkOperation(LARGE_NUMBER, "large number sqrt", sqrt_check);
    }
    {
        INFO("Проверка sin");
        checkOperation(SMALL_NUMBER, "small number sin", sin_check);
        checkOperation(LARGE_NUMBER, "large number sin", sin_check);
    }
    {
        INFO("Проверка arcsin");
        checkOperation(SMALL_NUMBER, "small number arcsin", asin_check);
        checkOperation(LARGE_NUMBER, "large number arcsin", asin_check);
    }
}