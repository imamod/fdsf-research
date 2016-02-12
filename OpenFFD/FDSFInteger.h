#pragma once

#include <vector>

namespace fdsfInteger{

    // Желаемая точность расчета
    const double epsilon = 1e-17;

    // Значение pi (Лучше брать из библиотеки mpfr)
    const double PI = 3.141592653589793238463;

    // Задает линейно-тригонометрическое распределение 
    void setLinTrigMesh(std::vector<double> &y0, std::vector<double> &x0, 
                        std::vector<double> &Y, std::vector<double> &X,
                        int N_base);

    // Вычисление Г-функции
    // TODO: сделать для полуцелых индексов
    // TODO: лучше подрубить boost
    double factorial(double k);

    // Вычисляет значение функции ФД индекса k = 1, 2, 3 в точке x при заданном
    // значении t
    double FermiDirakFunction(double t, double x, double k);

    // Вычисляет значение функции ФД индекса k = 1, 2, 3 в точке x по схеме 
    // Горнера при x <= -0.1. N - число членов в схеме Горнера для дрстижения 
    // машинной точности
    double Gorner(double x, int N, double k);

    // Вычисляет значение функции ФД индекса k = 1, 2, 3 в точке x по формуле 
    // Гаусса-Кристоффеля с 5 узлами на равномерной сетке. 
    // (сеточно Гауссов метод) N = 5
    // TODO: std::function
    double FDGK5(double(*Ft)(double, double, double), double x, double T, double k, int N);

    // Сгущение по Ричардсону результата работы функции FDGK5 
    double Richardson_mesh_refinement(double x, double t, double k, int N, int p);
} // namespace fdsfInteger