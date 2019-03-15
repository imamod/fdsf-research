#include "Common.h"
#include "FileSys.h"

namespace {

    const BmpReal X_LEFT = 0;

    void checkAccuracy(BmpReal k, BmpReal x, BmpReal left, BmpReal right) {
        setPreciseOutput();
        std::cout << "k = " << k << ", x = " << x << ", d = " << abs(left/right - 1) << std::endl;
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
