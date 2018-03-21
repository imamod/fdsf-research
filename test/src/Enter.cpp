#include "Common.h"
#include <iomanip>

namespace enter {
    // ������ ���������� �����
    size_t number(size_t upperBound) {
        std::cout << std::endl << "Enter �: ";
        size_t expected = 0;
        std::wcin >> expected;
        REQUIRE(expected < upperBound);
        return expected;
    }

    // ������� ������
    std::string string(const std::string& description) {
        std::cout << "Enter " << description << ":" << std::endl;
        std::string out;
        std::cin >> out;
        return out;
    }

    // ������� ������
    BmpReal index() {
        std::cout << "Select index:" << std::endl;
        for (size_t i = 0; i < ALL_INDICES.size(); ++i) {
            std::cout << "�" << i << " : " << std::setw(4) << ALL_INDICES[i] << std::endl;
        }
        size_t selected = enter::number(ALL_INDICES.size());
        return ALL_INDICES[selected];
    }

    // ������� �����
    BmpReal number(const std::string& description) {
        std::cout << "Enter " << description << " :" << std::endl;
        BmpReal x;
        std::cin >> x;
        return x;
    }
}
