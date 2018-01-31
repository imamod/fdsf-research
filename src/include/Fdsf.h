#pragma once

#include "BasicTypes.h"

namespace fdsf {

    // Значение Ik(0) для индекса k = 0,1,2,3
    const static BmpReal I_k_0[] = { log(BmpReal(2)),
                                      PI*PI / 12.0, 
                                      1.8030853547393952,// значение из статьи
                                      7.0*PI*PI*PI*PI / 120.0 
                                    };

    // Вычисление Г-функции
    BmpReal factorial(BmpReal k);

    // Вычисляет значение функции ФД индекса k = 1, 2, 3 в точке x при заданном
    // значении t. Функция представлена в виде t^k /(exp(x)+exp(t))
    BmpReal fermi_dirak_integer(BmpReal t, BmpReal x, BmpReal k);
#if 0
    BmpReal fermi_dirak_half_integer(BmpReal ksi, BmpReal x, BmpReal k);
#endif

    // Вычисляет значение функции ФД индекса k = 1, 2, 3 в точке x по схеме 
    // Горнера при x <= -0.1. N - число членов в схеме Горнера для достижения 
    // машинной точности
    BmpReal Gorner(BmpReal x, BmpReal k);

    // Алгоритм выбора T_max длякаждого x (а не фиксированные для всех)
    // NB: для новой формы записи подынтегральной функции
    BmpReal get_T_max(BmpReal X, int k);

    // Вычисляет значение функции ФД индекса k = 1, 2, 3 в точке x по формуле 
    // Гаусса-Кристоффеля с 5 узлами на равномерной сетке. 
    // По статье, для достижения требуемой точности на отрезке -0.1<x<0 нужно
    // брать для k = 0 : Tmax = 50, 
    //           k = 1 : Tmax = 60,
    //           k = 2 : Tmax = 75, 
    //           k = 3 : Tmax = 100. 
    // TODO: std::function
    BmpReal gauss_christoffel_method(BmpReal(*f)(BmpReal, BmpReal, BmpReal), 
                 BmpReal x, BmpReal T, BmpReal k, int N);

    // Сгущение по Ричардсону результата работы функции gauss_christoffel_method
    BmpReal richardson_method(BmpReal x, BmpReal k, BmpReal t = 0);

    struct integration_segment_values {
        size_t n; // Текущий отрезок интегрирования
        size_t N; // Общее число отрезков интегрирования
    };

    // for k = -3/2
    BmpReal fermi_dirak_m3half(BmpReal ksi, BmpReal x,
        BmpReal k, BmpReal a, integration_segment_values isv);

    // for others half-integer k
    BmpReal fermi_dirak_half_integer(BmpReal ksi, BmpReal x,
        BmpReal k, BmpReal a, integration_segment_values isv);

    // Formula Euler-Maclaurin
    BmpReal euler_maclaurin_method(BmpReal x, const BmpReal k, int N);

} // namespace fdsf
