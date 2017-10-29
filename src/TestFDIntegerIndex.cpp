#include "TestCommon.h"
#include "Fdsf.h"
#include <unordered_map>

namespace {
    // ��� �����
    BmpVector computeIntegral(const BmpVector& x, size_t k) {
        // std::cout << "      I1(0) : " << fdsf::PI*fdsf::PI / 12 << std::endl;
        /**
         * ������������ �������� t ��� ������� k, ��������� ����������������
         * bmp_real t = fdsf::get_T_max(X.at(i), k);
         */
        std::unordered_map<size_t, bmp_real> K_T_MAP = {
            { 1, 60 },
            { 2, 75 },
            { 3, 100 },
            { 4, 120 }
        };
        // �����, �� ������� ������� �� ����� �������
        bmp_real x_div = bmp_real(-0.1);
        BmpVector I;
        for (int i = 0; i < x.size(); i++) {
            auto value = x[i] > x_div ? fdsf::richardson_method(x[i], K_T_MAP[k], k)
                                         : fdsf::Gorner(x[i], k);
                //I.push_back(GornerSchemeForPrecesionY( x0[i], N));
            I.push_back(value);
            std::cout << "x0: " << x[i] << " I_base: " << I[i] << std::endl;
        }
        return I;
    }
}

TEST_CASE("k_1") {
    std::cout.precision(std::numeric_limits<bmp_real>::max_digits10);
    //int N_gorner = 260, k = 1;
    //int N_gorner = 214, k = 2;
    //int N_gorner = 165, k = 3;
    const size_t k = 1;
    // TODO: tests for N = 2 - 6
    const size_t N_base = 2;

    // ������ �������� ��������� � ������� �����
    BmpVector x0, X, y0, Y;
    fdsf::SetLinearTrigonometricGrid(y0, x0, Y, X, N_base);

    // ������ ����� ������� � �������� ��������� �� ��������� �����
    BmpVector I_base = computeIntegral(x0, k);
    BmpVector I_additional = computeIntegral(X, k);
    //GetValue_w(I_base, I_additional, y0, x0, Y, X, k);
}
