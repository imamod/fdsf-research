#include "TestCommon.h"

namespace {
    void inverse(const BmpMatrix& A) {
        CMatrix matrix(A);
        auto inversedMatrix = matrix.inverse();
        std::cout << "Printing inversed matrix: " << std::endl;
        matrix.print(inversedMatrix);
        std::cout << "--------------------------" << std::endl;
    }
}

TEST_CASE("Matrix_Tests") {
    std::cout << "Begining with matrixes" << std::endl;

    const BmpMatrix test1 = { { 3, 4 }, { 5, 7 } };
    const BmpMatrix test2 = { { 2, 5, 7 }, { 6, 3, 4 }, {5, -2, -3} };
    const BmpMatrix test3 = { { 2, 1, 0, 0 }, { 3, 2, 0, 0 }, { 1, 1, 3, 4 }, { 2, -1, 2, 3 } };

    const auto tests = { test1, test2, test3 };
    for (const auto& item : tests) {
        inverse(item);
    }

    BmpMatrix A, A_inv;
    // BmpVector B, x, delta_base(y0.size(), 0), delta_add(Y.size(), 0);
    BmpVector a, b;
    //object.fill_matrix(N_base, I_base, y0, B, A);

    //A_inv = inverse(A);
    //object.find_coefficients(A_inv, B, a, b, N_base);
    //printResultToFile(a, k, "a"); printResultToFile(b, k, "b");

    //GetApproxomateValues(a, b, y0, Y, I_additional, I_base, delta_base, delta_add, N_base);
    //printResultToFile(delta_base, k, "delta_base"); printResultToFile(delta_add, k, "delta_add");
}
