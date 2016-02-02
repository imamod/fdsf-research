#pragma once

#include <math.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>

// ���������� �-������� TODO: ������� ��� ��������� ��������
double factorial(double k);

// ������� ��
double func_Fermi(double t, double x, double k);

// ����� �������
double Gorner(double x, int N, double k);

// ����������� ����� � �������������� ������ ������-����������� (������� ������� �����) N = 3
double regular_mesh_with_three_points(double(*Ft)(double, double, double), double x, double T, double k, int N);

// ����������� ����� � �������������� ������ ������-����������� (������� ������� �����) N = 5
double regular_mesh_with_five_points(double(*Ft)(double, double, double), double x, double T, double k, int N);

// �������� �� ���������� �� �������-�������� ������
double Richardson_mesh_refinement(double x, double t, double k, int N, int p);

// ���������������� ����� � �������������� ������ ������-�����������
double method_quasiregular_mesh_with_GK(double(*Ft)(double, double, double), double x, double k, int N);

// ���������������� �����
double method_quasiregular_mesh(double(*Ft)(double, double, double), double x, double k, int N);