#pragma once

// ���������� �-������� TODO: ������� ��� ��������� ��������
double factorial(double k);

// ������� ��
double func_Fermi(double t, double x, double k);

// ����� �������
double Gorner(double x, int N, double k);

// ����������� ����� � �������������� ������ ������-����������� (������� ������� �����) N = 5
double regular_mesh_with_five_points(double(*Ft)(double, double, double), double x, double T, double k, int N);

// �������� �� ���������� �� �������-�������� ������
double Richardson_mesh_refinement(double x, double t, double k, int N, int p);
