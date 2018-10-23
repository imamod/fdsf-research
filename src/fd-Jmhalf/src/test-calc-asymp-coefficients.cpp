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

    // Вычисление Ck x-> +inf
    BmpVector calculateCRight(int N) {
        Logger log("calculateCRight");
        log.info(std::to_string(-pi()*pi() / 24));
        BmpReal k(-_1 / 2);
        BmpReal prod = 1, kPairsProd = (k + 1);
        BmpVector result = { 1 }, A = { 1 };//, -pi()*pi() / 24 };
        for (int q = 1; q <= N; ++q) {
            BmpReal mul = 2 * (_1 - pow(2, _1 - 2 * q));
            prod *= kPairsProd*(kPairsProd - 1);
            kPairsProd -= 2;
            A.push_back(mul*zetaFunction(2*q)*prod);
        }
       // print::vector(A);
        for (int n = 1; n <= N; ++n) {
            BmpReal sum = 0;
            for (int q = 0; q <= n; ++q) {
                sum += A.at(q)*A.at(n - q);
            }
            result.push_back(sum);
        }
        return result;
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
    INFO("Вычисление коэффициентов асимптотического ряда для интегральной ФД");
    BmpVector coeff = calculateCRight(12);
    // TODO: setPreciseOutput();
    filesys::writeFile("iffdRight.txt", coeff);
}