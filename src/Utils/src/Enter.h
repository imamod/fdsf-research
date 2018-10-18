#pragma once

namespace enter {
    // ¬вести корректный номер
    size_t number(size_t upperBound);
    // —читать число
    BmpReal number(const std::string& description);
    // —читать строку
    std::string string(const std::string& description);
    // —читать индекс
    BmpReal index();
}
