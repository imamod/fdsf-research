/*
 * Внешний подключаемый заголовочный файл
 */
#pragma once

#include "BasicTypes.h"

namespace fdsf {
    /* Функции ФД целого индекса */
    BmpReal fd_0(BmpReal x);
    BmpReal fd_1(BmpReal x);
    BmpReal fd_2(BmpReal x);
    BmpReal fd_3(BmpReal x);
    BmpReal fd_4(BmpReal x);

    /* Функции ФД полуцелого индекса */
    BmpReal fd_m3half(BmpReal x);
    BmpReal fd_m1half(BmpReal x);
    BmpReal fd_1half(BmpReal x);
    BmpReal fd_3half(BmpReal x);
    BmpReal fd_5half(BmpReal x);
    BmpReal fd_7half(BmpReal x);
}
