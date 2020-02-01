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

        // ”становить сетку базовых узлов по линейно-тригонометрическому закону
        virtual void setLinearTrigonometricGrid() {
            const BmpReal alpha = 2 / (2 + pi());
            const BmpReal one = BmpReal(1);
            const BmpReal num2 = BmpReal(2);
            const BmpReal x_star = BmpReal(0);
            const BmpReal y_star = BmpReal(log(1 + exp(x_star)));
            BmpReal baseSize = BmpReal(2 * baseNCount() - 1);
            // «адаютс€ базовые узлы интерпол€ции
            BmpVector baseGrid;
            for (size_t j = 1; j <= baseSize; j++) {
                baseGrid.push_back(y_star / num2*(num2 * alpha*j / baseSize
                    + (one - alpha)*(one - cos(pi()*j / baseSize))));
            }
            setBaseGrid(baseGrid);
            // «адаютс€ дополнительные точки
            setAdditionalDots();
        }

    };

    // TODO: модуль прецизионных аппроксимаций
    BmpReal preciseValuek12(BmpReal y) {
        BmpVector a = {
            1,
            0.4230390367951766,
            0.1150948861362855,
            0.0211596194006563,
            0.0028109730001318,
            0.0002646450259842,
            0.0000163590289177,
            0.0000005184504861,
        };

        BmpVector b = {
            1,
            0.1301458179817451,
            0.0444027175790622,
            0.0040809614661287,
            0.0005222693769483,
            0.0000264461085702,
            0.0000009662169431,
        };

        BmpReal nom = 0;
        for (auto n = 0; n < a.size(); ++n) {
            nom += a.at(n) * pow(y, n);
        }
        BmpReal denom = 0;
        for (auto n = 0; n < b.size(); ++n) {
            denom += b.at(n) * pow(y, n);
        }
        BmpReal z = nom / denom;
        BmpReal k = fdsf::index::P1_HALF;
        BmpReal result = factorial(k)*y*pow(z, k);
        return result;
    }

    BmpVector calcZ(const BmpVector& dots) {
        BmpVector z;
        for (auto y : dots) {
            auto x = log( exp(y) - 1);
            z.push_back(preciseValuek12(y) / fdsf::global_approx::low_temp::fd_1half(x));
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
