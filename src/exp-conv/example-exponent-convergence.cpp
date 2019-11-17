#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "BasicService.h"
#include "Fdsf-legacy.h"
#include "Constants.h"

#include <iostream>

namespace {
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
    BmpReal I_prec = fdsf::Gorner(x, k);
    std::cout << I_base << std::endl;
    std::cout << I_prec << std::endl;
}

TEST_CASE("ExpConverge") {
    std::cout.precision(std::numeric_limits<BmpReal>::max_digits10);
    // Для статьи о сверхстепенной сходимости
    // TODO:
    epc::Richardson(func_demo, 0, pi());
    /*epc::Richardson(func_cos, 0, 1);
    epc::Richardson(func_exp_cos, 0, 1);
    epc::Richardson(func_exp_sin, 0, 1);
    */
}