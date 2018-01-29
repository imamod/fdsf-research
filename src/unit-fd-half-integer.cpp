#include "TestCommon.h"
#include "Fdsf.h"
#ifdef HIGH_PRECISION
#include "Tables.h"
#endif

/**
 * TODO:
 *   * TEST_CASES : halfinteger, quadratures, approximations
 *   * unnamed namespace for static functions
 *   * problem with creating data, think about useful logic of folder creation with experimental data
 *   * common constants (to uppercase and to one place)
 */

BmpVector calculate_series_part(BmpReal k, BmpVector& X) {
    using namespace fdsf;
    BmpVector coeff_A = { pow(PI, 2) / 6.0,
        pow(PI, 4) / 90.0,
        pow(PI, 6) / 945.0,
        pow(PI, 8) / 9450.0,
        pow(PI, 10) / 93555.0,
        691.0 * pow(PI, 12) / 638512875.0
    };

    BmpVector series_value;

    for (size_t i = 0; i < X.size(); i++) {
        BmpReal coeff_C = 1;
        BmpReal nom = k + 1;
        BmpReal series_sum = BmpReal(1.0);
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

static BmpReal get_assympt_value(BmpReal x, BmpReal k) {
    BmpVector I_minus, I, series_part, X = { x };
    series_part = calculate_series_part(k, X);
    I_minus.push_back(fdsf::richardson_method(-x, k));
    I.push_back(I_minus[0] + series_part[0]);
    //return I[0];
    return series_part[0];
}

static BmpReal get_series_value(BmpReal x, BmpReal k) {
    BmpReal series_value = 0;
    auto N = log(fdsf::epsilon) / (x);

    for (size_t n = 1; n < N; ++n) {
        series_value += pow(-1.0, n - 1) * exp(n*x) / pow(n, k + 1);
        //std::cout << series_value << std::endl;
    }

    return fdsf::factorial(k)*series_value;
}

// *****************************************************************************
// Функции работы с прецизионными аппроксимациями
// *****************************************************************************
BmpVector computeIntegral(BmpVector x, BmpReal k) {
    BmpVector I;
    for (size_t i = 0; i < x.size(); ++i) {
        I.push_back(fdsf::richardson_method(x[i], k));
        //std::cout << "x0: " << x[i] << " I: " << I[i] << std::endl;
    }
    return I;
}

void calculateTableForK(BmpReal k, size_t N, const std::string& prefix = "test/", GRID_TYPE gridType = GRID_TYPE::LIN_TRIG_RIGHT) {
    // Расчет значения интеграла в базовых узлах
    Grid grid(N);
    //TODO: to GRID class
    switch (gridType) {
        case GRID_TYPE::LINEAR:
            grid.setLinearGrid();
            break;
        case GRID_TYPE::TRIGONOMETRIC:
            break;
        case GRID_TYPE::LIN_TRIG_LEFT:
            grid.setLinearTrigonometricGrid();
            break;
        case GRID_TYPE::LIN_TRIG_RIGHT:
            grid.setLinearTrigonometricGridRight();
            break;
        case GRID_TYPE::SHIFTED:
            break;
        default:
            break;
    }

    BmpVector y0 = grid.base();
    BmpVector Y = grid.additional();

    // Расчет интеграла
    BmpVector I_base = computeIntegral(grid.xByY(y0), k);
    BmpVector I_additional = computeIntegral(grid.xByY(Y), k);

    //std::string absFilename = filesys::createDirectory(12, it, "test/");
    //std::string absFilename = filesys::createDirectory(12, N, "test_x5_left/");
    //std::string absFilename = filesys::createDirectory(12, N, "test/");
    //std::string absFilename = filesys::createDirectory(12, it, "test_x5/");
    std::string absFilename = filesys::createDirectory(12, N, prefix);
    filesys::writeFile(absFilename + "y0.txt", y0);
    filesys::writeFile(absFilename + "I_base.txt", I_base);
    filesys::writeFile(absFilename + "I_add.txt", I_additional);
    filesys::writeFile(absFilename + "Y.txt", Y);
}

TEST_CASE("halfint") {
    SECTION("comp_kostya_and_precise") {
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

    SECTION("check_negative_quadrature_values") {
        std::cout.precision(std::numeric_limits<BmpReal>::max_digits10);
        const BmpReal k = BmpReal(1.0 / 2.0);
        BmpReal x = BmpReal(-0.1), I, I_precise;
        I = fdsf::richardson_method(x, k);
        I_precise = get_series_value(x, k);
        filesys::writeFile("../../test/test.txt", { I, I_precise });
        // TODO: добавить точность отдельно для double и mp
        REQUIRE(abs(I - I_precise) < 1e-17);
    }

    SECTION("calculate_asimpt_value") {
        BmpVector X, Y;
        const BmpReal k = BmpReal(1.0 / 2.0);
        BmpReal h = 0.1;
        BmpVector I, I_minus, series_part;

        //проверка на идиота при х = 30
        X.push_back(30.0);
        //X.push_back(log(exp(Y[0]) - 1));
        series_part = calculate_series_part(k, X);
        I_minus.push_back(fdsf::richardson_method(-X[0], k));
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
            I_minus.push_back(fdsf::richardson_method(-X[i], k));
            I.push_back(I_minus[i] + series_part[i]);
        }
#endif
        //printResultToFile(I, k, "Asimpt_check");
    }

    SECTION("shiftgrid") {
        const BmpReal k = BmpReal(1.0 / 2.0);
        //const std::vector<size_t> N_base = { 3, 5, 7, 9 };
        //const std::vector<size_t> N_base{ 3 };
        //const std::vector<size_t> N_base{ 5 };
        //const std::vector<size_t> N_base{ 7 };
        const std::vector<size_t> N_base{ 9 };
        // Расчет значения интеграла в базовых узлах
        for (const auto& it : N_base) {
            Grid grid(it);
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
            //const BmpReal TAU(1);
            //const BmpReal TAU(0.75);
            grid.shiftLinTrigGrid(DELTA_K12_N9, TAU);
            BmpVector y0 = grid.base();
            BmpVector Y = grid.additional();
            BmpVector x0 = grid.xByY(y0);
            BmpVector X = grid.xByY(Y);
#ifdef HIGH_PRECISION
            BmpVector I_base = computeIntegral(x0, k);
            BmpVector I_additional(I_k12_n3);
#else
            // Расчет интеграла
            BmpVector I_base = computeIntegral(x0, k);
            BmpVector I_additional = computeIntegral(X, k);
#endif

            std::string absFilename = filesys::createDirectory(12, it, "testShift/");
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
    }

    SECTION("SolveSystem") {
        const BmpReal k = BmpReal(1.0 / 2.0);
        const std::vector<size_t> N_base{ 9 };
        for (const auto& it : N_base) {
            Grid grid(it);
            grid.setLinearTrigonometricGridRight();
            BmpVector y0 = grid.base();
            BmpVector Y = grid.additional();
            BmpVector x0 = grid.xByY(y0);
            BmpVector X = grid.xByY(Y);

            // Расчет интеграла
            BmpVector I_base = computeIntegral(x0, k);
            filesys::writeFile("I_base.txt", I_base);
            //BmpVector I_additional = computeIntegral(X, k);
            //filesys::writeFile("I_add.txt", I_additional);
            //BmpVector I_base = filesys::readFile("I_base.txt");
            //test::printVector(I_base, true);
            BmpVector I_additional = filesys::readFile("I_add.txt");
#if 0
            BmpVector E = solveRightApproximationSystem(k, it, y0, I_base);
            // Раскладываем вектор Е в коэффициенты a b аппроксимации
            BmpVector a(E.begin(), E.begin() + it);
            BmpVector b(E.begin() + it, E.end());

            // Вывод результата
            std::cout << "-----a-----" << std::endl;
            test::printVector(a, true);
            std::cout << "-----b-----" << std::endl;
            test::printVector(b, true);
            filesys::writeFile("a.txt", a);
            filesys::writeFile("b.txt", b);
            getchar();
#endif
        }
    }
}

TEST_CASE("quadrature") {
    SECTION("m32") {
        const BmpReal k = BmpReal(-3.0 / 2.0);
        const std::vector<size_t> N_base = { 3, 5, 7, 9 };
        // Расчет значения интеграла в базовых узлах
        for (const auto& it : N_base) {
            calculateTableForK(k, it);
        }
    }

    SECTION("m12") {
        const BmpReal k = BmpReal(-1.0 / 2.0);
        const std::vector<size_t> N_base = { 3, 5, 7, 9 };
        // Расчет значения интеграла в базовых узлах
        for (const auto& it : N_base) {
            calculateTableForK(k, it);
        }
    }

    SECTION("12") {
        const BmpReal k = BmpReal(1.0 / 2.0);
        const std::vector<size_t> N_base = { 3, 5, 7, 9 };
        // Расчет значения интеграла в базовых узлах
        for (const auto& it : N_base) {
            calculateTableForK(k, it);
        }
    }

    SECTION("32") {
        const BmpReal k = BmpReal(3.0 / 2.0);
        const std::vector<size_t> N_base = { 3, 5, 7, 9 };
        // Расчет значения интеграла в базовых узлах
        for (const auto& it : N_base) {
            calculateTableForK(k, it);
        }
    }

    SECTION("52") {
        const BmpReal k = BmpReal(5.0 / 2.0);
        const std::vector<size_t> N_base = { 3, 5, 7, 9 };
        // Расчет значения интеграла в базовых узлах
        for (const auto& it : N_base) {
            calculateTableForK(k, it);
        }
    }

    SECTION("72") {
        const BmpReal k = BmpReal(7.0 / 2.0);
        const std::vector<size_t> N_base = { 3, 5, 7, 9 };
        // Расчет значения интеграла в базовых узлах
        for (const auto& it : N_base) {
            calculateTableForK(k, it);
        }
    }
}

/*******************************************************************************
 * Проверялась чебышевская сетка, пока не удалять
 *******************************************************************************/
void chebyshevBaseNodes(BmpVector &y_base,
                        BmpVector &x_base,
                        BmpVector &Y,
                        BmpVector &X, int N_base) {
    using namespace fdsf;
    const size_t n_additional = 11;
    const BmpReal one = BmpReal(1);
    const BmpReal num2 = BmpReal(2); //if integer
    const BmpReal x_star = BmpReal(3);
    const BmpReal y_star = BmpReal(log(1 + exp(x_star))); // if half-integer
    BmpReal baseSize = BmpReal(N_base); // if poly approximation

    const BmpReal y_star_inv = 1 / (y_star * y_star);

    // Задаются базовые узлы интерполяции
    for (size_t j = 1; j <= baseSize; j++) {
        y_base.push_back(y_star_inv*pow(sin(PI*j / (2 * baseSize)), 2));
        std::cout << y_base[j - 1] << " ";
    }

    // Задаются дополнительные точки
    Y.push_back(y_base[0] / n_additional);

    for (size_t i = 1; i < n_additional; i++) {
        Y.push_back(Y[i - 1] + y_base[0] / n_additional);
    }

    for (size_t index = 1; index < y_base.size(); index++) {
        for (size_t i = 0; i < n_additional; i++) {
            Y.push_back(Y.back() + (y_base[index] - y_base[index - 1]) / n_additional);
        }
    }

    // Разворачиваем y
    std::reverse(y_base.begin(), y_base.end());
    for (size_t j = 0; j < baseSize; j++) {
        y_base[j] = 1.0 / pow(y_base[j], 0.5);
        x_base.push_back(log(exp(y_base[j]) - one));
    }

    std::reverse(Y.begin(), Y.end());
    //std::cout << Y.size() << std::endl;
    for (size_t j = 0; j < Y.size(); j++) {
        Y[j] = 1.0 / pow(Y[j], 0.5);
        X.push_back(log(exp(Y[j]) - one));
    }
}


TEST_CASE("test_chebyshev_base_nodes") {
    BmpReal k = 0.5;
    BmpVector x0, X, y0, Y;
    const size_t N_base = 1;
    // Расчет значения интеграла в базовых узлах
    chebyshevBaseNodes(y0, x0, Y, X, N_base);
    // Расчет интеграла
    BmpVector I_base = computeIntegral(x0, k);
    BmpVector I_additional = computeIntegral(X, k);

    filesys::writeFile("I_base", I_base);
    filesys::writeFile("I_add", I_additional);
    filesys::writeFile("y0", y0);
    filesys::writeFile("Y", Y);
}

TEST_CASE("calculate_all_k") {
    BmpVector k = { -0.5, 0.5, 1.5, 2.5, 3.5 };
    BmpReal left_start = -5.0, right_start = 3;
    BmpReal left_end = 5.0, right_end = 30;
    BmpReal span = 0.1;
    BmpVector x_left, x_right;
    for (auto i = left_start; i < left_end; i += span) {
        x_left.push_back(i);
    }
    for (auto i = right_start; i < right_end; i += span) {
        x_right.push_back(i);
    }
    //std::cout << x_left.size() << std::endl;
    for (auto index : k) {
        BmpVector I_left;
        for (auto item : x_left) {
            I_left.push_back(fdsf::richardson_method(item, index));
            //filesys::writeFile("I", I_left);
        }
        BmpVector I_right;
        for (auto item : x_right) {
            I_right.push_back(fdsf::richardson_method(item, index));
            //filesys::writeFile("I_right", I_right);
        }
        //std::cout << I.size() << std::endl;
    }

   /* BmpVector I_m32;
    for (auto item : x_left) {
        I_m32.push_back(fdsf::fd_3mhalf(item));
    }
    filesys::writeFile("I_m32", I_m32);

    BmpVector I_m32_right;
    for (auto item : x_right) {
        I_m32_right.push_back(fdsf::fd_3mhalf(item));
    }
    filesys::writeFile("I_right_m32", I_m32_right);*/
}
