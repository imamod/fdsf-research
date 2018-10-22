#include "Common.h"
#include "FileSys.h"

/*******************************************************************************
 * Проверялась чебышевская сетка, пока не удалять
 *******************************************************************************/
namespace {
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
            y_base.push_back(y_star_inv*pow(sin(pi()*j / (2 * baseSize)), 2));
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
}


TEST_CASE("test_chebyshev_base_nodes") {
    BmpReal k = 0.5;
    BmpVector x0, X, y0, Y;
    const size_t N_base = 1;
    // Расчет значения интеграла в базовых узлах
    chebyshevBaseNodes(y0, x0, Y, X, N_base);
    // Расчет интеграла
    BmpVector I_base = compute::halfInteger(x0, k);
    BmpVector I_additional = compute::halfInteger(X, k);

    filesys::writeFile("I_base", I_base);
    filesys::writeFile("I_add", I_additional);
    filesys::writeFile("y0", y0);
    filesys::writeFile("Y", Y);
}