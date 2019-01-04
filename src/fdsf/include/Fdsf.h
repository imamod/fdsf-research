/*
 * Внешний подключаемый заголовочный файл
 */
#pragma once

#include "BasicTypes.h"

namespace fdsf {
    /* Поддерживаемые индексы функций ФД */
    namespace index {
        const BmpReal M3_HALF = -3.0 / 2;
        const BmpReal M1_HALF = -1.0 / 2;
        const BmpReal ZERO = 0;
        const BmpReal P1_HALF = 1.0 / 2;
        const BmpReal P1 = 1;
        const BmpReal P3_HALF = 3.0 / 2;
        const BmpReal P2 = 2;
        const BmpReal P5_HALF = 5.0 / 2;
        const BmpReal P3 = 3;
        const BmpReal P7_HALF = 7.0 / 2;
        const BmpReal P4 = 4;
    }

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
