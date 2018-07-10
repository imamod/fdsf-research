#pragma once

#include "BasicService.h"

// Значение pi
inline const BmpReal pi() {
    return boost::math::constants::pi<BmpReal>();
}

// Константа Эйлера
inline const BmpReal eilerConst() {
    return 0.5772156649015325;
}

// Константа j
inline const BmpReal j() {
    return 0.76740941382814898;
}

// Кузьмина дисер чек
inline const BmpReal jFromKuzmina() {
    return 0.95024075;
}
