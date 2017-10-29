#include "TestCommon.h"
#include "Fdsf.h"

std::vector<bmp_real> calculate_series_part(bmp_real k, std::vector<bmp_real>& X) {
    using namespace fdsf;
    std::vector<bmp_real> coeff_A = { pow(PI, 2) / 6.0,
        pow(PI, 4) / 90.0,
        pow(PI, 6) / 945.0,
        pow(PI, 8) / 9450.0,
        pow(PI, 10) / 93555.0,
        691.0 * pow(PI, 12) / 638512875.0
    };

    std::vector<bmp_real> series_value;

    for (size_t i = 0; i < X.size(); i++) {
        bmp_real coeff_C = 1;
        bmp_real nom = k + 1;
        bmp_real series_sum = bmp_real(1.0);
        for (size_t j = 0; j < coeff_A.size(); j++) {
            coeff_C *= nom*(nom - 1); // По асимптотической формуле парное добавление множителей, поэтому далее отнимаем 2
            series_sum += 2.0 * (1.0 - pow(2.0, 1.0 - 2 * (j + 1))) * pow(X[i], (-2.0)*(j + 1))*coeff_A[j] * coeff_C;
            //std::cout << "A(j) = " << coeff_A[j] << ": series_sum = " << series_sum << std::endl;
            std::cout << "C(" << j + 1 << ") = " << 2.0 * (1.0 - pow(2.0, 1.0 - 2 * (j + 1))) * coeff_A[j] * coeff_C << std::endl;
            nom -= 2;
        }
        series_sum *= pow(X[i], k + 1) / (k + 1);
        series_value.push_back(series_sum);
    }

    std::cout << PI*PI*PI*PI*PI*PI / 945.0 << std::endl;

    return series_value;
}

static bmp_real get_assympt_value(bmp_real x, bmp_real k) {
    std::vector<bmp_real> I_minus, I, series_part, X = { x };
    bmp_real t = 0, a = 0;
    series_part = calculate_series_part(k, X);
    I_minus.push_back(fdsf::richardson_method(-x, t, k, a));
    I.push_back(I_minus[0] + series_part[0]);
    //return I[0];
    return series_part[0];
}

static bmp_real get_series_value(bmp_real x, bmp_real k) {
    bmp_real series_value = 0;
    auto N = log(fdsf::epsilon) / (x);

    for (size_t n = 1; n < N; ++n) {
        series_value += pow(-1.0, n - 1) * exp(n*x) / pow(n, k + 1);
        //std::cout << series_value << std::endl;
    }

    return fdsf::factorial(k)*series_value;
}

TEST_CASE("comp_kostya_and_precise") {
    const bmp_real k = bmp_real(1.0 / 2.0);
    bmp_real x = bmp_real(-0.1), I, I_kostya, I_precise;
    bmp_real t = 0, a = 0;
#if 0
    I = fdsf::richardson_method(x, t, k, a);
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
    I = fdsf::richardson_method(x, t, k, a);
    //I_kostya = fdsf::fd_half(x);
    I_precise = get_assympt_value(x, k);
    std::cout << "I_quadrature: " << I << std::endl;
    //std::cout << "I_kostya: " << I_kostya << std::endl;
    std::cout << "I_precise: " << I_precise << std::endl;
    std::cout << "delta = " << I / I_precise - 1 << std::endl;
}

TEST_CASE("check_negative_quadrature_values") {
    std::cout.precision(std::numeric_limits<bmp_real>::max_digits10);
    const bmp_real k = bmp_real(1.0 / 2.0);
    bmp_real x = bmp_real(-0.1), I, I_precise;
    bmp_real t = 0, a = 0;
    I = fdsf::richardson_method(x, t, k, a);
    I_precise = get_series_value(x, k);
    std::cout << "I_quadrature: " << I << std::endl;
    std::cout << "I_precise: " << I_precise << std::endl;
    REQUIRE(abs(I - I_precise) < 1e-17);
}

TEST_CASE("calculate_asimpt_value") {
    std::vector<bmp_real> X, Y;
    const bmp_real k = bmp_real(1.0 / 2.0);
    bmp_real t = 0, a = 0, h = 0.1;
    std::vector<bmp_real> I, I_minus, series_part;

    //проверка на идиота при х = 30
    X.push_back(30.0);
    //X.push_back(log(exp(Y[0]) - 1));
    series_part = calculate_series_part(k, X);
    I_minus.push_back(fdsf::richardson_method(-X[0], t, k, a));
    I.push_back(I_minus[0] + series_part[0]);
    I.push_back(series_part[0]);
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

    series_part = calculate_series_part(k, X);
    for (size_t i = 0; i < X.size(); i++) {
        I_minus.push_back(fdsf::richardson_method(-X[i], t, k, a));
        I.push_back(I_minus[i] + series_part[i]);
    }
#endif
    //printResultToFile(I, k, "Asimpt_check");
}