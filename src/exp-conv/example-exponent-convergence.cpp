#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "BasicService.h"
#include "Gamma.h"
#include "Constants.h"

#include <iostream>

namespace fdsf {
// Сгущение по Ричардсону по сеточно-Гауссову методу
BmpReal richardson_method(BmpReal x, BmpReal k, BmpReal t = 0);
}

namespace {

    // Получение необходимого числа членов для схемы Горнера произвольной 
    // заданной точности
    BmpReal get_N_for_Gorner(BmpReal x, BmpReal k) {
        return log(fdsf::epsilon) / x;
    }

    // Вычисляет значение функции ФД индекса k = 1, 2, 3 в точке x по схеме 
    // Горнера при x <= -0.1. N - число членов в схеме Горнера для достижения 
    // машинной точности
    BmpReal Gorner(BmpReal x, BmpReal k) {
        //size_t N = ceil(get_N_for_Gorner(x, k));
        size_t N = (size_t)(get_N_for_Gorner(x, k)) + 1;
        BmpReal exp_x = exp(x);
        BmpReal sum = 1.0 / pow(N, k + 1);

        for (size_t i = N - 1; i > 0; i--) {
            sum = 1 / pow(i, k + 1) - exp_x * sum;
        }

        return sum * factorial(k) * exp(x);
    }

    BmpReal func_cos(BmpReal x) {
        return 1.0 / (2 - cos(x));
    }

    BmpReal func_demo(BmpReal x) {
        size_t n = 0;
        BmpReal a = 1.75;
        return (a*a - 1)*pow(a, n)*cos(n*x) / pow(1 - 2 * a*cos(x) + a*a, 1);
    }

    BmpReal func_exp_sin(BmpReal x) {
        return exp(sin(x));
    }

    BmpReal func_exp_cos(BmpReal x) {
        return exp(cos(x));
    }
}

// Сравнение схемы Горнера и метода трапеций epc
TEST_CASE("GornerVsTrapz") {
    //TODO: setPreciseOutput();
    BmpReal k = 1.0 / 2;
    BmpReal x = -1;
    BmpReal I_base = fdsf::richardson_method(x, k);
    BmpReal I_prec = Gorner(x, k);
    std::cout << I_base << std::endl;
    std::cout << I_prec << std::endl;
}

TEST_CASE("ExpConverge") {
    std::cout.precision(std::numeric_limits<BmpReal>::max_digits10);
    // Для статьи о сверхстепенной сходимости
    // TODO:
    for (BmpReal p : {0, 1, 2}) {
        for (BmpReal c : { exp(0.3162), exp(1), exp(3.162) }) {
            auto func_demo_check = [c, p](BmpReal x) {
                return (c*c - 1)*pow(c, p)*cos(p*x) / pow(1 - 2 * c*cos(x) + c*c, 1);
            };
            epc::Richardson(func_demo_check, 0, pi(), c, p);
        }
    }
    // epc::Richardson(func_demo, 0, pi());
    /*epc::Richardson(func_cos, 0, 1);
    epc::Richardson(func_exp_cos, 0, 1);
    epc::Richardson(func_exp_sin, 0, 1);
    */
}