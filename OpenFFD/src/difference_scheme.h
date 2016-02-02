#pragma once

#include <math.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>

// Вычисление Г-функции TODO: сделать для полуцелых индексов
double factorial(double k);

// Функция ФД
double func_Fermi(double t, double x, double k);

// Схема Горнера
double Gorner(double x, int N, double k);

// Равномерная сетка с использованием формул Гаусса-Кристофелля (сеточно Гауссов метод) N = 3
double regular_mesh_with_three_points(double(*Ft)(double, double, double), double x, double T, double k, int N);

// Равномерная сетка с использованием формул Гаусса-Кристофелля (сеточно Гауссов метод) N = 5
double regular_mesh_with_five_points(double(*Ft)(double, double, double), double x, double T, double k, int N);

// Сгущение по Ричардсону по сеточно-Гауссову методу
double Richardson_mesh_refinement(double x, double t, double k, int N, int p);

// Квазиравномерная сетка с использованием формул Гаусса-Кристофелля
double method_quasiregular_mesh_with_GK(double(*Ft)(double, double, double), double x, double k, int N);

// Квазиравномерная сетка
double method_quasiregular_mesh(double(*Ft)(double, double, double), double x, double k, int N);