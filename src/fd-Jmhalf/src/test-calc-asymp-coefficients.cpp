#include "Common.h"
#include "SeriesLogNDivSqrN.h"
#include "FileSys.h"
#include "ZetaFunction.h"
#include "Logger.h"

namespace {

    // Вычисление константы j
    inline const BmpReal calculateConstj() {
        SeriesLogNDivSqrN example;
        BmpReal value = (1 - 2 * log(2) / 3 - eilerConst() / 3);
        BmpReal jConst = pow(pi(), 2)*value / 2 - example.get();
        return jConst;
    }

    void piPowerCheck(const BmpReal& value) {
        BmpReal piPower = 0.5;
        BmpReal upperBound = enter::number(UINT32_MAX);
        // TODO: setPreciseOutput();
        while (piPower < upperBound) {
            std::cout << " power = " << piPower << ": " << pow(pi(), piPower) / value << std::endl;
            piPower += 0.5;
        }
    }

    /**
     * Вычисление Ck x-> +inf для [I-0.5]^2, для вычисления коэффициентов
     * интегральной ф-ии каждый поделить на (n-1)
     */
    BmpVector calculateIm2Coeff(int N) {
        Logger log("calculateCRight");
        BmpReal k(-_1 / 2);
        BmpReal prod = 1, kPairsProd = (k + 1);
        BmpVector result = { 1 }, A = { 1 };//, -pi()*pi() / 24 };
        for (int q = 1; q <= N; ++q) {
            BmpReal mul = 2 * (_1 - pow(2, _1 - 2 * q));
            prod *= kPairsProd*(kPairsProd - 1);
            kPairsProd -= 2;
            A.push_back(mul*zetaFunction(2*q)*prod);
        }
        //print::vector(A);
        for (int n = 1; n <= N; ++n) {
            BmpReal sum = 0;
            for (int q = 0; q <= n; ++q) {
                sum += A.at(q)*A.at(n - q);
            }
            result.push_back(sum);
        }
        return result;
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
            value += C_N.at(n) / ((n - 1) * pow(x, 2 * n));
        }
        return value;
    }

    const BmpReal CONST_J = 0.76740941382814898;

    // Обратный ход вычисления j
    BmpReal asymptotic(BmpReal x) {
        BmpReal seriesValue = seriesSum(x);
        BmpReal bracket = 1 + (CONST_J - pi()*pi()*log(x) / 6) / (x*x) - seriesValue;
        BmpReal value = 2 * x*x* bracket;
        return value;
    }


    // Вычисляет коэффициенты интегральной ф-ии по коэффициентам Im2half
    BmpVector Jm2Coeff() {
        const BmpVector IM2_HALF = {
            1.00000000000000000,
            -0.82246703342411309,
            -3.38226010534730559,
            -56.74866767632004638,
            -2076.43981697169328982,
            -133516.62391908300924115,
            -13363920.49546855688095093,
            -1924202279.42978835105895996,
            -376996608458.57202148437500000,
            -96469021655492.73437500000000000,
            -31243036135798104.00000000000000000,
            -12492545181655248896.00000000000000000,
            -6044381261816933646336.00000000000000000,
        };
        BmpVector result{ 1, -pi()*pi()/6 };
        for (size_t n = 2; n < IM2_HALF.size(); ++n) {
            result.push_back(-IM2_HALF.at(n)/(n-1));
        }
        return result;
    }

    void checkAccuracy(BmpReal x, BmpReal left, BmpReal right) {
        setPreciseOutput();
        std::cout << " x = " << x << ", d = " << left / right - 1 << std::endl;
    }
}

TEST_CASE("seriesCheck") {
    SECTION("convergence") {
        std::cout << "N=100: " << SeriesLogNDivSqrN(100).get() << std::endl;
        std::cout << "N=200: " << SeriesLogNDivSqrN(200).get() << std::endl;
        std::cout << "N=300: " << SeriesLogNDivSqrN().get() << std::endl;
        std::cout << "N=10^5 " << SeriesLogNDivSqrN(100000).get() << std::endl;
    }
    SECTION("pi-power") {
        SeriesLogNDivSqrN series;
        piPowerCheck(series.get());
    }
    SECTION("check-eiler-plus-sum") {
        SeriesLogNDivSqrN series;
        // TODO: setPreciseOutput();
        BmpReal eilerPlusSeriesSum = pow(pi(), 2)*eilerConst() / 6 + series.get();
        std::cout << eilerPlusSeriesSum << std::endl;
        piPowerCheck(eilerPlusSeriesSum);
    }
}

TEST_CASE("calculate") {
    {
        INFO("Вычисление коэффициентов асимптотического ряда для [Im2]^2 ФД");
        nlohmann::json coeff = calculateIm2Coeff(12);
        // TODO: setPreciseOutput();
        filesys::writeFile("Im2.json", coeff);
    }
    {
        INFO("Вычисление коэффициентов асимптотического ряда для интегральной ФД");
        nlohmann::json coeff = Jm2Coeff();
        filesys::writeFile("Jm2.json", coeff);
    }
}

TEST_CASE("accuracy") {
    SECTION("quad_asym") {
        const BmpVector QUAD_30_60 = {
            1790.3529161501417,
            2439.8437531948134,
            3189.4031422037642,
            4039.0147557204718,
            4988.6674939821187,
            6038.353463537884,
            7188.0668485663755,
        };
        size_t i = 0;
        for (auto x : { 30, 35, 40, 45, 50, 55, 60 }) {
           checkAccuracy(x, asymptotic(x), QUAD_30_60.at(i++));
        }
    }
}

TEST_CASE("asympt_borders") {
    BmpVector coeff = {
        1.00000000000000000,
        - 0.82246703342411309,
        - 3.38226010534730559,
        - 56.74866767632004638,
        - 2076.43981697169328982,
        - 133516.62391908300924115,
        - 13363920.49546855688095093,
        - 1924202279.42978835105895996,
        - 376996608458.57202148437500000,
        - 96469021655492.73437500000000000,
        - 31243036135798104.00000000000000000,
        - 12492545181655248896.00000000000000000,
        - 6044381261816933646336.00000000000000000,
    };
    setPreciseOutput();
    BmpVector firstEq;
    for (int i = coeff.size() - 1; i > 0; --i) {
        BmpReal relation = coeff.at(i) / coeff.at(i - 1);
        //std::cout << "C(" << i << ")/C(" << i-1 << ") = " << relation << std::endl;
        BmpReal value = 2 * sqrt(relation*(i-1)/(i-2));
        firstEq.push_back(value);
        std::cout << " n = " << i - 1 << " : x = " << value << std::endl;
    }
    BmpVector secondEq;
    for (int i = coeff.size() - 1; i > 0; --i) {
        BmpVector precisions;
        std::cout << " n = " << i - 1 << " : x = ";
        for (auto eps : { 1e-8, 1e-16, 1e-19 }) {
            BmpReal value = pow(4*(i - 1)*abs(coeff.at(i)) / (3.0 * eps), 1.0 / (2*i));
            std::cout << value << " ";
        }
        std::cout << std::endl;
    }
}