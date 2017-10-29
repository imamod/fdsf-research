#pragma once

#include "BasicTypes.h"

namespace fdsf {

    // Значение Ik(0) для индекса k = 0,1,2,3
    const static bmp_real I_k_0[] = { log(bmp_real(2)),
                                      PI*PI / 12.0, 
                                      1.8030853547393952,// значение из статьи
                                      7.0*PI*PI*PI*PI / 120.0 
                                    };

    // Задает линейно-тригонометрическое сетку в базовых узлах, плюс 10(?)
    // дополнительных точек между каждой парой базовых узлов.
    void SetLinearTrigonometricGrid(std::vector<bmp_real> &y_base,
                                    std::vector<bmp_real> &x_base,
                                    std::vector<bmp_real> &Y,
                                    std::vector<bmp_real> &X, size_t N_base);

    // Вычисление Г-функции
    bmp_real factorial(bmp_real k);

    // Вычисляет значение функции ФД индекса k = 1, 2, 3 в точке x при заданном
    // значении t. Функция представлена в виде t^k /(exp(x)+exp(t))
    bmp_real fermi_dirak_integer(bmp_real t, bmp_real x, bmp_real k);
#if 0
    bmp_real fermi_dirak_half_integer(bmp_real ksi, bmp_real x, bmp_real k);
#endif

    // Вычисляет значение функции ФД индекса k = 1, 2, 3 в точке x по схеме 
    // Горнера при x <= -0.1. N - число членов в схеме Горнера для достижения 
    // машинной точности
    bmp_real Gorner(bmp_real x, bmp_real k);

    // Алгоритм выбора T_max длякаждого x (а не фиксированные для всех)
    // NB: для новой формы записи подынтегральной функции
    bmp_real get_T_max(bmp_real X, int k);

    // Вычисляет значение функции ФД индекса k = 1, 2, 3 в точке x по формуле 
    // Гаусса-Кристоффеля с 5 узлами на равномерной сетке. 
    // По статье, для достижения требуемой точности на отрезке -0.1<x<0 нужно
    // брать для k = 0 : Tmax = 50, 
    //           k = 1 : Tmax = 60,
    //           k = 2 : Tmax = 75, 
    //           k = 3 : Tmax = 100. 
    // TODO: std::function
    bmp_real gauss_christoffel_method(bmp_real(*f)(bmp_real, bmp_real, bmp_real), 
                 bmp_real x, bmp_real T, bmp_real k, int N);

    // Сгущение по Ричардсону результата работы функции gauss_christoffel_method
    bmp_real richardson_method(bmp_real x, bmp_real t, bmp_real k, bmp_real a = 0);

    struct integration_segment_values {
        size_t n; // Текущий отрезок интегрирования
        size_t N; // Общее число отрезков интегрирования
    };

    // for k = -3/2
    bmp_real fermi_dirak_m3half(bmp_real ksi, bmp_real x,
        bmp_real k, bmp_real a, integration_segment_values isv);

    // for others half-integer k
    bmp_real fermi_dirak_half_integer(bmp_real ksi, bmp_real x,
        bmp_real k, bmp_real a, integration_segment_values isv);

    // Formula Euler-Maclaurin 
    bmp_real euler_maclaurin_method(bmp_real x, const bmp_real k, int N, bmp_real& a);

} // namespace fdsf
