#include "Common.h"
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


// *****************************************************************************
// Функции работы с прецизионными аппроксимациями
// *****************************************************************************

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
    BmpVector I_base = compute::halfInteger(grid.xByY(y0), k);
    BmpVector I_additional = compute::halfInteger(grid.xByY(Y), k);

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
            BmpVector I_base = compute::halfInteger(x0, k);
            filesys::writeFile("I_base.txt", I_base);
            //BmpVector I_additional = compute::halfInteger(X, k);
            //filesys::writeFile("I_add.txt", I_additional);
            //BmpVector I_base = filesys::readFile("I_base.txt");
            BmpVector I_additional = filesys::readFile("I_add.txt");
#if 0
            BmpVector E = solveRightApproximationSystem(k, it, y0, I_base);
            // Раскладываем вектор Е в коэффициенты a b аппроксимации
            BmpVector a(E.begin(), E.begin() + it);
            BmpVector b(E.begin() + it, E.end());

            // Вывод результата
            std::cout << "-----a-----" << std::endl;
            print::vector(a, true);
            std::cout << "-----b-----" << std::endl;
            print::vector(b, true);
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
