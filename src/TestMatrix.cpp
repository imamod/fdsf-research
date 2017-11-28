#include "TestCommon.h"

TEST_CASE("MatrixTests") {
    // map { data: expected }
    const std::map<BmpMatrix, BmpMatrix> TESTS = {
        { { { 3, 4 }, { 5, 7 } },
          { { 7, -4 }, { -5, 3 } }
        },
        { { { 2, 5, 7 }, { 6, 3, 4 }, { 5, -2, -3 } },
          { { 1, -1, 1 }, { -38, 41, -34 }, { 27, -29, 24 } }
        },
        { { { 2, 1, 0, 0 }, { 3, 2, 0, 0 }, { 1, 1, 3, 4 }, { 2, -1, 2, 3 } },
          { { 2, -1, 0, 0 }, { -3, 2, 0, 0 }, { 31, -19, 3, -4 }, { -23, 14, -2, 3 } }
        },
    };

    for (const auto& item : TESTS) {
        auto inversedMatrix = CMatrix(item.first).inverse();
        //REQUIRE(item.second == inversedMatrix);
        test::printMatrix(inversedMatrix);
    }
    // TODO: move to integer tests
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

TEST_CASE("ExampleGrid") {
    size_t N = 5;
    Grid grid(N);
    grid.setLinearGrid();
    auto base = grid.base();
    auto additional = grid.additional();
    std::cout << "Base: ";
    test::printVector(base);
    std::cout << "Additional: ";
    test::printVector(additional);
}
