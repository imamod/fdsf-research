#include "Common.h"
#include <iomanip>

namespace enter {
    // ¬вести корректный номер
    size_t number(size_t upperBound) {
        std::cout << std::endl << "Enter є: ";
        size_t expected = 0;
        std::wcin >> expected;
        REQUIRE(expected < upperBound);
        return expected;
    }

    // —читать строку
    std::string string(const std::string& description) {
        std::cout << "Enter " << description << ":" << std::endl;
        std::string out;
        std::cin >> out;
        return out;
    }

    // —читать индекс
    BmpReal index() {
        std::cout << "Select index:" << std::endl;
        for (size_t i = 0; i < ALL_INDICES.size(); ++i) {
            std::cout << "є" << i << " : " << std::setw(4) << ALL_INDICES[i] << std::endl;
        }
        size_t selected = enter::number(ALL_INDICES.size());
        return ALL_INDICES[selected];
    }

    // —читать число
    BmpReal number(const std::string& description) {
        std::cout << "Enter " << description << " :" << std::endl;
        BmpReal x;
        std::cin >> x;
        return x;
    }
}
