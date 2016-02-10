#pragma once

#include <vector>

namespace fdsfInteger{

    // Желаемая точность расчета
    const double epsilon = 1e-17;

    // Значение pi (Лучше брать из библиотеки mpfr)
    const double PI = 3.141592653589793238463;

    // Смешанная сетка
    void setLinTrigMesh(std::vector<double> &y0, std::vector<double> &x0, std::vector<double> &Y, std::vector<double> &X, int N_base);

    // Вычисление Г-функции
    // TODO: сделать для полуцелых индексов
    // TODO: лучше подрубить boost
    double factorial(double k);

    // Функция ФД
    double func_Fermi(double t, double x, double k);

    // Схема Горнера
    double Gorner(double x, int N, double k);

    // Равномерная сетка с использованием формул Гаусса-Кристофелля (сеточно Гауссов метод) N = 5
    double regular_mesh_with_five_points(double(*Ft)(double, double, double), double x, double T, double k, int N);

    // Сгущение по Ричардсону по сеточно-Гауссову методу
    double Richardson_mesh_refinement(double x, double t, double k, int N, int p);
}; // namespace fdsfInteger