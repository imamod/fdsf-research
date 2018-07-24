#include "Common.h"

namespace print {

    // Распечатать вектор
    void vector(const BmpVector& v, bool printInColumn) {
        std::cout << std::setprecision(std::numeric_limits<BmpReal>::digits10 + 2) << std::endl;
        auto escapeSymbol = printInColumn ? "\n" : " ";
        for (auto const& it : v) {
            std::cout << it << escapeSymbol;
        }
        std::cout << std::endl;
    }

    // Распечатать матрицу
    void matrix(const BmpMatrix& m) {
        std::cout << std::setprecision(std::numeric_limits<BmpReal>::digits10 + 2) << std::endl;
        for (auto it : m) {
            print::vector(it);
        }
        std::cout << std::endl;
    }

    // Распечатать мапу (удобно для пар x, f(x))
    void map(const std::map<BmpReal, BmpReal>& m) {
        std::cout << std::setprecision(std::numeric_limits<BmpReal>::digits10 + 2) << std::endl;
        for (auto it : m) {
            std::cout << it.first << std::setw(4) << " : " << it.second << std::endl;
        }
        std::cout << std::endl;
    }
}

/* Установить расширенный потоковый вывод на 16 знаков */
void setPreciseOutput() {
    std::cout.precision(std::numeric_limits<BmpReal>::max_digits10);
}
