#include "Common.h"
#include "Filesys.h"

namespace {
    const BmpVector DELTA_K12_N3 = { 9.438492e-06, 6.092994e-06, 5.484654e-06,
        3.481478e-06, 6.460857e-07, 2.988557e-08 };
    const BmpVector DELTA_K12_N5 = { 3.815749e-09, 2.141089e-09, 2.207771e-09,
        3.366675e-09, 6.286565e-09, 1.155590e-08,
        1.381581e-08, 5.660865e-09, 4.789431e-10,
        4.265144e-12
    };
    const BmpVector DELTA_K12_N7 = { 5.360046e-12, 2.032152e-12, 1.448508e-12,
        1.552092e-12, 2.300937e-12, 4.358069e-12,
        1.000855e-11, 2.630696e-11, 7.034207e-11,
        1.614167e-10, 2.371363e-10, 1.542799e-10,
        3.149425e-11, 2.248868e-12 };
    const BmpVector DELTA_K12_N9 = { 3.996803e-15, 1.554312e-15, 9.992007e-16,
        5.551115e-16, 6.661338e-16, 6.661338e-16,
        1.443290e-15, 3.552714e-15, 7.327472e-15,
        2.431388e-14, 4.773959e-14, 3.599343e-13,
        1.252554e-12, 3.502754e-12, 5.019984e-12,
        2.974954e-12, 5.715428e-13, 3.841372e-14 };

    BmpVector deltaByN(size_t N) {
        switch (N) {
            case 3:
                return DELTA_K12_N3;
            case 5:
                return DELTA_K12_N5;
            case 7:
                return DELTA_K12_N7;
            case 9:
                return DELTA_K12_N9;
        }
    }
}

TEST_CASE("shiftgrid") {
    const BmpReal k = BmpReal(1.0 / 2.0);
    const std::vector<size_t> N_base_possible = { 3, 5, 7, 9 };
    std::cout << "Posible N_base: " << std::endl;
    for (size_t i = 0; i < N_base_possible.size(); ++i) {
        std::cout << "№" << i << " : " << N_base_possible[i] << std::endl;
    }
    const size_t N_base = enter::number(N_base_possible.size());
    // Расчет значения интеграла в базовых узлах
    Grid grid(N_base);
    grid.setLinearTrigonometricGridRight();
    // BmpVector delta = ;
     /* для N = 9 лин-триг сетка, не удалять !!!!!*/
     /*BmpVector delta = { 3.996803e-015,
                         1.554312e-015,
                         9.992007e-016,
                         6.661338e-016,
                         5.551115e-016,
                         1.110223e-015,
                         1.554312e-015,
                         2.442491e-015,
                         2.886580e-015,
                         1.088019e-014,
                         1.136868e-013,
                         2.143841e-013,
                         1.398881e-012,
                         3.310685e-012,
                         5.213607e-012,
                         2.875700e-012,
                         6.026291e-013,
                         4.263256e-014
                     };*/
                     // TAU = 0.75
    BmpVector delta = { 3.996803e-15, 1.554312e-15, 9.992007e-16, 5.551115e-16,
        6.661338e-16, 6.661338e-16, 1.443290e-15, 3.552714e-15,
        7.327472e-15, 2.431388e-14, 4.773959e-14, 3.599343e-13,
        1.252554e-12, 3.502754e-12, 5.019984e-12, 2.974954e-12,
        5.715428e-13, 3.841372e-14,
    };
    const BmpReal TAU(0.5);
    // N = 3
    BmpReal tau = enter::number(3);
    //const BmpReal TAU(1);
    //const BmpReal TAU(0.75);
    grid.shiftLinTrigGrid(deltaByN(N_base), TAU);
    BmpVector y0 = grid.base();
    BmpVector Y = grid.additional();
    BmpVector x0 = grid.xByY(y0);
    BmpVector X = grid.xByY(Y);
#ifdef HIGH_PRECISION
    BmpVector I_base = compute::halfInteger(x0, k);
    BmpVector I_additional(I_k12_n3);
#else
    // Расчет интеграла
    BmpVector I_base = compute::halfInteger(x0, k);
    BmpVector I_additional = compute::halfInteger(X, k);
#endif

    std::string absFilename = filesys::createDirectory(12, N_base, "testShift/");
    filesys::writeFile(absFilename + "y0.txt", y0);
    filesys::writeFile(absFilename + "I_base.txt", I_base);
    filesys::writeFile(absFilename + "I_add.txt", I_additional);
    filesys::writeFile(absFilename + "Y.txt", Y);
    //getchar();
#if 0
    BmpVector E = solveRightApproximationSystem(k, it, y0, I_base);
    // Раскладываем вектор Е в коэффициенты a b аппроксимации
    BmpVector a(E.begin(), E.begin() + it);
    BmpVector b(E.begin() + it, E.end());
    filesys::writeFile(absFilename + "a.txt", a);
    filesys::writeFile(absFilename + "b.txt", b);

    BmpVector I_app = approximateValueRight(a, b, y0, k);
    filesys::writeFile(absFilename + "I_app.txt", I_app);
    BmpVector I_app_add = approximateValueRight(a, b, Y, k);
    filesys::writeFile(absFilename + "I_app_add.txt", I_app_add);
#endif
}