#pragma once

#include <vector>

namespace fdsfInteger{

    // Желаемая точность расчета
    const double epsilon = 1e-17;

    // Значение pi (Лучше брать из библиотеки mpfr)
    const double PI = 3.141592653589793238463;

    // Значение Ik(0) для индекса k = 0,1,2,3
    const static double I_k_0[] = { log(2),
                                    PI*PI / 12.0, 
                                    1.8030853547393952,// значение из статьи
                                    7.0*PI*PI*PI*PI / 120.0 };

    // Задает линейно-тригонометрическое сетку в базовых узлах, плюс 10(?)
    // дополнительных точек между каждой парой базовых узлов.
    void SetLinesrTrigonometricGrid(std::vector<double> &y_base, 
                                    std::vector<double> &x_base, 
                                    std::vector<double> &Y, 
                                    std::vector<double> &X, int N_base);

    // Вычисление Г-функции
    // TODO: сделать для полуцелых индексов
    // TODO: лучше подрубить boost
    double factorial(double k);

    // Вычисляет значение функции ФД индекса k = 1, 2, 3 в точке x при заданном
    // значении t. Функция представлена в виде t^k /(exp(x)+exp(t))
    double FermiDirakFunction(double t, double x, double k);

    // Вычисляет значение функции ФД индекса k = 1, 2, 3 в точке x по схеме 
    // Горнера при x <= -0.1. N - число членов в схеме Горнера для достижения 
    // машинной точности
    double Gorner(double x, int N, double k);

    // Вычисляет значение функции ФД индекса k = 1, 2, 3 в точке x по формуле 
    // Гаусса-Кристоффеля с 5 узлами на равномерной сетке. 
    // (сеточно Гауссов метод) N = 5
    // По статье, для достижения требуемой точности на отрезке -0.1<x<0 нужно
    // брать для k = 0 : Tmax = 50, 
    //           k = 1 : Tmax = 60,
    //           k = 2 : Tmax = 75, 
    //           k = 3 : Tmax = 100.
    // TODO: Перенести алгоритм автоматического выбора Tmax.
    // TODO: std::function
    double FDGK5(double(*Ft)(double, double, double), 
                 double x, double T, double k, int N);

    // Сгущение по Ричардсону результата работы функции FDGK5
    double Richardson_mesh_refinement(double x, double t, double k, int N);

    // Вычисляет значение функции ФД индекса k=1 в точке x
    double FD_I1(double x);

    // Вычисляет значение функции ФД индекса k=2 в точке x
    double FD_I2(double x);

    // Вычисляет значение функции ФД индекса k=3 в точке x
    double FD_I3(double x);
} // namespace fdsfInteger