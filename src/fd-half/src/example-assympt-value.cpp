#include "Common.h"
#include "Gamma.h"
#include "Constants.h"
#include <iostream>

namespace {

    // Cтруктура хранения x_min и N_max (см. препринт2, табл 4)
    struct Limits {
        BmpReal x_min;
        int N_max;
    };

    // Получение x_min и N_max для конкретного индекса
    Limits limits(BmpReal k) {
        if (k == -1.5) {
            return { 44, 11 };
        } else if (k == -0.5) {
            return { 39, 10 };
        } else if (k == 0.5) {
            return { 35, 10 };
        } else if (k == 1.5) {
            return{ 33, 10 };
        } else if (k == 2.5) {
            return { 30, 10 };
        } else if (k == 3.5) {
            return { 29, 10 };
        }
        throw std::invalid_argument("Unsuppported index");
    }

    // TODO сделать тест по вычислению коэффициентов в отдельном test-файле
    // Перенести этот кодв AsymptoticSeries класс удалить
    // Получает коэффициенты асимптотического ряда для конкретного k
    BmpVector coefficents(BmpReal k) {
        BmpVector A;
        // TODO: сделать гибко, через limits
        //constexpr int N = 10;
        constexpr int N = 12;

        for (int n = 1; n <= N; ++n) {
            const BmpReal constMember = 2.0 - pow(2, 2 - 2*n);
            BmpReal prod = 1;
            for (int p = 1; p <= 2 * n; ++p) {
                prod *= k + 2 - p;
            }
            BmpReal result = constMember*dzetaFunction(2 * n)*prod;
            A.push_back(result);
        }
        return A;
    }

    // Вычислить значение ФД для x >= x_min
    // Схема Горнера для асимптотического ряда (см. формулу (32) препринт 2)
    BmpReal calculate(BmpReal k, BmpReal x) {
        // Получаем коэффициенты для конкретного k
        BmpVector A = coefficents(k);
        // Получаем предельные значения для каждого k
        Limits data = limits(k);
        BmpReal x_2_m1 = 1.0 / (x*x);
        BmpReal sum = A.back() * x_2_m1;
        // TODO: когда поймем, какое N_max, использовать его
        for (int n = A.size() - 2; n > -1; --n) {
            sum = (A.at(n) + sum) * x_2_m1;
            if (n > 3) {
                std::cout << "A(" << n+2 << ")/A(" << n+1 << ") = " << (A.at(n + 1)/A.at(n))*x_2_m1 << std::endl;
            }
        }
        // Главный член асимптотики
        BmpReal mainPart = pow(x, k + 1) / (k + 1);
        return mainPart * (sum + 1);
    }
}

TEST_CASE("calculate") {
    setPreciseOutput();
    SECTION("m12") {
        INFO("Вычисление значения функции ФД индекса k = -1/2 при x >= x_min");
        BmpReal I_x_min = calculate(-0.5, 39);
        std::cout << "k = -0.5 I(39) = " << I_x_min << std::endl;
    }
    SECTION("12") {
        INFO("Вычисление значения функции ФД индекса k = 1/2 при x >= x_min");
        BmpReal I_x_min = calculate(0.5, 35);
        std::cout << "k = 0.5 I(35) = " << I_x_min << std::endl;
    }
    SECTION("32") {
        INFO("Вычисление значения функции ФД индекса k = 3/2 при x >= x_min");
        BmpReal I_x_min = calculate(1.5, 33);
        std::cout << "k = 1.5 I(33) = " << I_x_min << std::endl;
    }
    SECTION("52") {
        INFO("Вычисление значения функции ФД индекса k = 5/2 при x >= x_min");
        BmpReal I_x_min = calculate(2.5, 30);
        std::cout << "k = 2.5 I(30) = " << I_x_min << std::endl;
    }
    SECTION("72") {
        INFO("Вычисление значения функции ФД индекса k = 7/2 при x >= x_min");
        BmpReal I_x_min = calculate(3.5, 29);
        std::cout << "k = 3.5 I(29) = " << I_x_min << std::endl;
    }
}

TEST_CASE("dzeta") {
    for (size_t n = 2; n <= 24; n += 2) {
        std::cout << n / 2 << " "<< dzetaFunction(n) << std::endl;
    }
}