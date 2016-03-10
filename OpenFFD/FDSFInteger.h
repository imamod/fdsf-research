#pragma once

#include <vector>
#include "mpreal.h"

namespace fdsf{

    // Желаемая точность расчета
    const int digits = 25;
    const mpfr::mpreal epsilon = pow(10,-digits);

    // Значение pi (Лучше брать из библиотеки mpfr)
    //const mpfr::mpreal PI = 3.141592653589793238463;
    const mpfr::mpreal PI = mpfr::const_pi();

    // Значение Ik(0) для индекса k = 0,1,2,3
    const static mpfr::mpreal I_k_0[] = { log(2),
                                    PI*PI / 12.0, 
                                    1.8030853547393952,// значение из статьи
                                    7.0*PI*PI*PI*PI / 120.0 };

    // Задает линейно-тригонометрическое сетку в базовых узлах, плюс 10(?)
    // дополнительных точек между каждой парой базовых узлов.
    void SetLinearTrigonometricGrid(std::vector<mpfr::mpreal> &y_base, 
                                    std::vector<mpfr::mpreal> &x_base, 
                                    std::vector<mpfr::mpreal> &Y, 
                                    std::vector<mpfr::mpreal> &X, int N_base);

    // Вычисление Г-функции
    // TODO: сделать для полуцелых индексов
    // TODO: лучше подрубить boost
    mpfr::mpreal factorial(mpfr::mpreal k);

    // Вычисляет значение функции ФД индекса k = 1, 2, 3 в точке x при заданном
    // значении t. Функция представлена в виде t^k /(exp(x)+exp(t))
    mpfr::mpreal FermiDirakFunction(mpfr::mpreal t, mpfr::mpreal x, mpfr::mpreal k);

    // Вычисляет значение функции ФД индекса k = 1, 2, 3 в точке x по схеме 
    // Горнера при x <= -0.1. N - число членов в схеме Горнера для достижения 
    // машинной точности
    mpfr::mpreal Gorner(mpfr::mpreal x, int N, mpfr::mpreal k);

    // Алгоритм выбора T_max длякаждого x (а не фиксированные для всех)
    // NB: для новой формы записи подынтегральной функции
    mpfr::mpreal get_T_max(mpfr::mpreal X, double k);

    // Вычисляет значение функции ФД индекса k = 1, 2, 3 в точке x по формуле 
    // Гаусса-Кристоффеля с 5 узлами на равномерной сетке. 
    // По статье, для достижения требуемой точности на отрезке -0.1<x<0 нужно
    // брать для k = 0 : Tmax = 50, 
    //           k = 1 : Tmax = 60,
    //           k = 2 : Tmax = 75, 
    //           k = 3 : Tmax = 100. 
    // TODO: std::function
    mpfr::mpreal FDGK5(mpfr::mpreal(*Ft)(mpfr::mpreal, mpfr::mpreal, mpfr::mpreal), 
                       mpfr::mpreal x, mpfr::mpreal T, mpfr::mpreal k, int N);

    // Сгущение по Ричардсону результата работы функции FDGK5
    mpfr::mpreal Richardson_mesh_refinement(mpfr::mpreal x, mpfr::mpreal t, mpfr::mpreal k, int N);
    
    namespace integer {
        // Вычисляет значение функции ФД индекса k=1 в точке x
        mpfr::mpreal FD_I1(mpfr::mpreal x);

        // Вычисляет значение функции ФД индекса k=2 в точке x
        mpfr::mpreal FD_I2(mpfr::mpreal x);

        // Вычисляет значение функции ФД индекса k=3 в точке x
        mpfr::mpreal FD_I3(mpfr::mpreal x);
    } // namespace integer

} // namespace fdsf

mpfr::mpreal GornerSchemeForPrecesionY(mpfr::mpreal x, int N);