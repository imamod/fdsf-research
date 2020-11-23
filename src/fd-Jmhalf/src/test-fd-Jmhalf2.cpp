#include "Common.h"
#include "Logger.h"

namespace {

    const BmpReal PI2_DIV_6 = pi()*pi()/6;

    const inline BmpReal fdJmhalf2X30() {
        return 1790.3529161501788;
    }

    const BmpVector C_N = {
        1.00000000000000000,
        -0.82246703342411309,
        -3.38226010534730559,
        -56.7486676763200464,
        -2076.43981697169329,
        -133516.623919083009,
        -13363920.4954685569,
        -1924202279.42978835,
        -376996608458.572022,
        -96469021655492.7344,
        -31243036135798104.0,
        -12492545181655248896.0,
        -6044381261816933646336.0,
    };

    // Вычисление суммы ряда в асимптотическом разложении
    BmpReal seriesSum(BmpReal x) {
        Logger logger("seriesSum");
        int N = 12;
        BmpReal value = 0;
        for (int n = N; n > 1; --n) {
            logger.info("n = " + std::to_string(n));
            value += C_N.at(n) / ((n - 1) * pow(x, 2*n));
        }
        return 2 * x * x * value;
    }

    // Обратный ход вычисления j
    BmpReal calculateConstJ(BmpReal x, BmpReal J_x) {
        Logger logger("calculateConstJ");
        BmpReal seriesValue = seriesSum(x);
        logger.info("seriesSum " + std::to_string(seriesValue));
        BmpReal value = 3 * (J_x - 2 * x*x + seriesValue) / (pi()*pi()) + log(x);
        return value;
    }
}

TEST_CASE("jRestore") {
    setPreciseOutput();
    const BmpReal x = 30.0;
    BmpReal delta = _1 - fdJmhalf2X30() / (2*x*x); // 0.005359...
    BmpReal jRestored = log(x) - 6*x*x*delta / (pi()*pi());
    std::cout << x << ") : " << jRestored << std::endl;
}

TEST_CASE("jCheck") {
    setPreciseOutput();
    SECTION("new") {
        const BmpVector QUAD_30_50 = {
            1790.3529161501417,
            2439.8437531948134,
            3189.4031422037642,
            4039.0147557204718,
            4988.6674939821187,
        };
        size_t i = 0;
        for (auto x : {30, 35, 40, 45, 50}) {
            std::cout << calculateConstJ(x, QUAD_30_50.at(i++)) << std::endl;
        }
    }
    SECTION("gsl") {
        const BmpVector QUAD_30_50 = {
            1790.352916150179,
            2439.843753194812,
            3189.4031422038,
            4039.014755720504,
            4988.667493982043,
        };
        size_t i = 0;
        for (auto x : { 30, 35, 40, 45, 50 }) {
            std::cout << calculateConstJ(x, QUAD_30_50.at(i++)) << std::endl;
        }
    }
}

TEST_CASE("check") {
    setPreciseOutput();
    std::cout << 0.76740941382814898 / PI2_DIV_6 << std::endl;
}

TEST_CASE("coeff") {
    setPreciseOutput();
    for (auto x : { pi()*pi() / 6, 5*pow(pi(), 4) / 144, 17*pow(pi(), 6) / 576}) {
        std::cout << x << std::endl;
    }
}