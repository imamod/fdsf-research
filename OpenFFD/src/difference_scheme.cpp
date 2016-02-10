// Implementation of difference_scheme.h
//
#include "matrix_helper.h"
#include "difference_scheme.h"
#include <math.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>

// Функция ФД
double func_Fermi(double t, double x, double k) 
{
    // testfunction for positive X
    //double u = pow(t,k)/(1 + exp(t-x));
    // Для расчета интегралом, Для схемы Горнера
    double u = pow(t, k) / (factorial(k)*(exp(x) + exp(t)));
    return u;
}

// Схема Горнера
double Gorner(double x, int N, double k)
{
    double exp_x = exp(x);
    double sum = 1.0 / pow(N, k + 1);

    for (int i = N - 1; i > 0; i--){
        sum = 1 / pow(i, k + 1) - exp_x*sum;
    }

    return sum;
}

// Вычисление Г-функции TODO: сделать для полуцелых индексов
double factorial(double k)
{
    if (k < 0)
        return 0;
    if (k == 0)
        return 1;
    else
        return k*factorial(k - 1);
}

// Равномерная сетка с использованием формул Гаусса-Кристофелля (сеточно Гауссов метод) n = 5 пятиточечная схема
double regular_mesh_with_five_points(double(*Ft)(double, double, double), double x, double T, double k, int N)
{
    double gamma_1_5 = (322.0 - 13.0*sqrt(70)) / 1800.0;
    double gamma_2_4 = (322.0 + 13.0*sqrt(70)) / 1800.0;
    double gamma_3 = 64.0 / 225.0;
    matrix_type::_vector t(5);

    double U = 0;

    for (int n = N-1; n >=0 ; n--) {
        t.at(0) = T*(n + 0.5 - 0.5*sqrt((35 + 2 * sqrt(70)) / 63)) / N;
        t.at(1) = T*(n + 0.5 - 0.5*sqrt((35 - 2 * sqrt(70)) / 63)) / N;
        t.at(2) = T*(n + 0.5) / N;
        t.at(3) = T*(n + 0.5 + 0.5*sqrt((35 - 2 * sqrt(70)) / 63)) / N;
        t.at(4) = T*(n + 0.5 + 0.5*sqrt((35 + 2 * sqrt(70)) / 63)) / N;
        U = U + T*(Ft(t.at(2), x, k)*gamma_3 + gamma_1_5*((Ft(t.at(0), x, k)) + Ft(t.at(4), x, k)) + gamma_2_4*((Ft(t.at(1), x, k)) + Ft(t.at(3), x, k)));
    }

    U = U / N;

    return U;
}

// Сгущение по Ричардсону по сеточно-Гауссову методу
double Richardson_mesh_refinement(double x, double t, double k, int N, int p)
{
    double current_accuracy;
    double I_n, I_2n, I;

    I_n = regular_mesh_with_five_points(&func_Fermi, x, t, k, N);
    do {
        I_2n = regular_mesh_with_five_points(&func_Fermi, x, t, k, 2 * N);
        current_accuracy = (I_2n - I_n) / (pow(2.0, p) - 1);
        I = I_2n + current_accuracy;
        I_n = I_2n;
        N = 2 * N;

    } while (abs(current_accuracy) > epsilon); // Фактическая точность 10^-16

    return I;
}