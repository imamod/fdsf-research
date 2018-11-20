#pragma once

#include "BasicTypes.h"

// Некоторые числовые константы
const BmpReal _1 = BmpReal(1);
const BmpReal _2 = BmpReal(2);
const BmpReal _5 = BmpReal(5);
const BmpReal _7 = BmpReal(7);
const BmpReal _691 = BmpReal(691);
const BmpReal _3617 = BmpReal(3617);
const BmpReal _43867 = BmpReal(43867);
const BmpReal _77683 = BmpReal(77683);
const BmpReal _155366 = BmpReal(155366);
const BmpReal _174611 = BmpReal(174611);
const BmpReal _854513 = BmpReal(854513);
const BmpReal _236364091 = BmpReal(236364091);

// Значение pi
inline const BmpReal pi() {
    return 4 * atan(1);
    //return boost::math::constants::pi<BmpReal>();
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
