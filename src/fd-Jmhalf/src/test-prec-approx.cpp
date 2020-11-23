#include "Common.h"
#include "json.hpp"
#include "Fdsf.h"


namespace {

    class GridJmhalf : public Grid {
        public:
            GridJmhalf(size_t N_base, size_t addNCount)
                : Grid(N_base, addNCount) {}

            // Установить сетку базовых узлов по линейно-тригонометрическому закону
            void setLinearTrigonometricGrid() {
                const BmpReal alpha = 2 / (2 + pi());
                const BmpReal one = BmpReal(1);
                const BmpReal num2 = BmpReal(2);
                const BmpReal x_star = BmpReal(5);
                const BmpReal y_star = BmpReal(log(1 + exp(x_star)));
                BmpReal baseSize = BmpReal(2 * baseNCount());
                // Задаются базовые узлы интерполяции
                for (size_t j = 1; j <= baseSize; j++) {
                    baseGrid().push_back(y_star / num2*(num2 * alpha*j / baseSize
                        + (one - alpha)*(one - cos(pi()*j / baseSize))));
                }

                // Задаются дополнительные точки
                setAdditionalDots();
            }

    };

}

TEST_CASE("left") {
    // Расчет базовых узлов
    GridJmhalf grid(4, 11);
    grid.setLinearTrigonometricGrid();
    fd;
}
