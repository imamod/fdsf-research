#include "Common.h"
#include "Gamma.h"
#include "AsymptoticSeries.h"
#include <iostream>

namespace {

    BmpReal get_assympt_value(BmpReal x, BmpReal k) {
        BmpVector I_minus, I, X = { x };
        AsymptoticSeries series(k, x);
        std::cout << "new " << series.get() << std::endl;
        I_minus.push_back(fdsf::richardson_method(-x, k));
        I.push_back(I_minus[0] + series.get());
        //return I[0];
        return series.get();
    }

    BmpReal get_series_value(BmpReal x, BmpReal k) {
        BmpReal series_value = 0;
        auto N = log(fdsf::epsilon) / (x);

        for (size_t n = 1; n < N; ++n) {
            series_value += pow(-1.0, n - 1) * exp(n*x) / pow(n, k + 1);
            //std::cout << series_value << std::endl;
        }

        return factorial(k)*series_value;
    }
}

TEST_CASE("calculate_asimpt_value") {
    BmpVector X, Y;
    const BmpReal k = BmpReal(1.0 / 2.0);
    BmpReal h = 0.1;
    BmpVector I, I_minus;

    //проверка на идиота при х = 30
    X.push_back(30.0);
    //X.push_back(log(exp(Y[0]) - 1));
    AsymptoticSeries series(k, X.front());
    I_minus.push_back(fdsf::richardson_method(-X[0], k));
    I.push_back(I_minus[0] + series.get());
    I.push_back(series.get());
#if 0
    Y.push_back(3.0);
    X.push_back(log(exp(Y[0]) - 1));
    size_t i = 1;
    while (true)
    {
        Y.push_back(Y[0] + i*h);
        X.push_back(log(exp(Y[i]) - 1));

        if (Y[i] > 30.0) {
            break;
        }

        i++;
    }

    for (size_t i = 0; i < X.size(); i++) {
        AsymptoticSeries series(k, X.at(i));
        I_minus.push_back(fdsf::richardson_method(-X[i], k));
        I.push_back(I_minus[i] + series.get(i));
    }
#endif
    //printResultToFile(I, k, "Asimpt_check");
}

TEST_CASE("check_negative_quadrature_values") {
    setPreciseOutput();
    const BmpReal k = BmpReal(1.0 / 2.0);
    BmpReal x = BmpReal(-0.1), I, I_precise;
    I = fdsf::richardson_method(x, k);
    I_precise = get_series_value(x, k);
    filesys::writeFile("../../test/test.txt", { I, I_precise });
    // TODO: добавить точность отдельно для double и mp
    REQUIRE(abs(I - I_precise) < 1e-17);
}

TEST_CASE("comp_kostya_and_precise") {
        const BmpReal k = BmpReal(1.0 / 2.0);
        BmpReal x = BmpReal(-0.1), I, I_kostya, I_precise;
#if 0
        I = fdsf::richardson_method(x, k);
        I_kostya = fdsf::fd_half(x);
        I_precise = get_series_value(x, k);
        std::cout << "x = -0.1" << std::endl;
        std::cout << "I_quadrature: " << I << std::endl;
        std::cout << "I_kostya: " << I_kostya << std::endl;
        std::cout << "I_precise: " << I_precise << std::endl;
        std::cout << "delta = " << I / I_precise - 1 << std::endl;
#endif

        x = 30.0;// +10E-8;
        std::cout << "x = " << x << std::endl;
        I = fdsf::richardson_method(x, k);
        //I_kostya = fdsf::fd_half(x);
        I_precise = get_assympt_value(x, k);
        std::cout << "I_quadrature: " << I << std::endl;
        //std::cout << "I_kostya: " << I_kostya << std::endl;
        std::cout << "I_precise: " << I_precise << std::endl;
        std::cout << "delta = " << I / I_precise - 1 << std::endl;
}
