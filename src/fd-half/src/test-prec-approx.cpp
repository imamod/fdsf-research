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

    };

    BmpVector calcZ(const BmpVector& dots) {
        BmpVector z;
        for (auto y : dots) {
            auto x = log( exp(y) - 1);
            z.push_back(fdsf::prec_approx::fd_1half(x) / fdsf::global_approx::low_temp::fd_1half(x));
        }
        return z;
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
