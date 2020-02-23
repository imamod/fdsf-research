#include "Common.h"
#include "Grid.h"
#include "Fdsf.h"
#include "FdIndex.h"
#include "Gamma.h"
#include "Filesys.h"

namespace {

    class GridIhalf : public Grid {
    public:
        GridIhalf(size_t N_base, size_t addNCount)
            : Grid(N_base, addNCount) {}

        // Установить сетку базовых узлов по линейно-тригонометрическому закону
        virtual void setLinearTrigonometricGrid() {
            const BmpReal alpha = 2 / (2 + pi());
            const BmpReal one = BmpReal(1);
            const BmpReal num2 = BmpReal(2);
            const BmpReal x_star = BmpReal(0);
            const BmpReal y_star = BmpReal(log(1 + exp(x_star)));
            BmpReal baseSize = BmpReal(2 * baseNCount() - 1);
            // Задаются базовые узлы интерполяции
            BmpVector baseGrid;
            for (size_t j = 1; j <= baseSize; j++) {
                baseGrid.push_back(y_star / num2*(num2 * alpha*j / baseSize
                    + (one - alpha)*(one - cos(pi()*j / baseSize))));
            }
            setBaseGrid(baseGrid);
            // Задаются дополнительные точки
            setAdditionalDots();
        }

       /* virtual void setLinearTrigonometricGridRight() {
            additional().clear();
            base().clear();
            const BmpReal alpha = 2 / (2 + pi());
            const BmpReal one = BmpReal(1);
            const BmpReal num2 = BmpReal(2);
            const BmpReal x_star = BmpReal(7);
            const BmpReal y_star = BmpReal(log(1 + exp(x_star)));
            BmpReal baseSize = BmpReal(2 * baseNCount());

            const BmpReal P(2);
            const BmpReal y_star_new = 1 / pow(y_star, P);
            // Задаются базовые узлы интерполяции
            BmpVector baseGrid;
            for (size_t j = 1; j <= baseSize; j++) {
                BmpReal linearPart = num2 * (one - alpha)*j / baseSize;
                BmpReal trigonometricPart = alpha*(one - cos(pi()*j / baseSize));
                baseGrid.push_back(pow(y_star_new / num2*(linearPart + trigonometricPart), 1.0 / P));
            }
            setBaseGrid(baseGrid);
            // Задаются дополнительные точки
            setAdditionalDots();
        }*/
    };

    BmpVector calcZ(const BmpVector& dots) {
        BmpVector z;
        for (auto y : dots) {
            auto x = log( exp(y) - 1);
            z.push_back(fdsf::prec_approx::fd_1half(x) / fdsf::global_approx::low_temp::fd_1half(x));
        }
        return z;
    }

    BmpVector calcIRight(const BmpVector& dots) {
        BmpVector I;
        for (auto y : dots) {
            auto x = log(exp(y) - 1);
            I.push_back(fdsf::fd_1half(x));
        }
        return I;
    }
}

TEST_CASE("left") {
    GridIhalf grid(4, 20);
    grid.setLinearTrigonometricGrid();
    setPreciseOutput();
    filesys::writeFile("y_base", grid.base());
    filesys::writeFile("y_add", grid.additional());
    filesys::writeFile("z_base", calcZ(grid.base()));
    filesys::writeFile("z_add", calcZ(grid.additional()));
}

TEST_CASE("right") {
    GridIhalf grid(5, 11);
    grid.setLinearTrigonometricGridRight();
    setPreciseOutput();
    filesys::writeFile("y_base_right", grid.base());
    filesys::writeFile("y_add_right", grid.additional());
    filesys::writeFile("I_base_right", calcIRight(grid.base()));
    filesys::writeFile("I_add_right", calcIRight(grid.additional()));
}

TEST_CASE("check") {
    setPreciseOutput();
    std::cout << fdsf::fd_m1half(0) << std::endl;
    std::cout << pow(fdsf::fd_m1half(0), 4) << std::endl;
}
