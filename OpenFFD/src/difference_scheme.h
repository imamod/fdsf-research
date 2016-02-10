#pragma once

// Вычисление Г-функции TODO: сделать для полуцелых индексов
double factorial(double k);

// Функция ФД
double func_Fermi(double t, double x, double k);

// Схема Горнера
double Gorner(double x, int N, double k);

// Равномерная сетка с использованием формул Гаусса-Кристофелля (сеточно Гауссов метод) N = 5
double regular_mesh_with_five_points(double(*Ft)(double, double, double), double x, double T, double k, int N);

// Сгущение по Ричардсону по сеточно-Гауссову методу
double Richardson_mesh_refinement(double x, double t, double k, int N, int p);
