#pragma once

#include <vector>

namespace fdsfInteger{

    // �������� �������� �������
    const double epsilon = 1e-17;

    // �������� pi (����� ����� �� ���������� mpfr)
    const double PI = 3.141592653589793238463;

    // ������ �������-������������������ ������������� 
    void setLinTrigMesh(std::vector<double> &y0, std::vector<double> &x0, 
                        std::vector<double> &Y, std::vector<double> &X,
                        int N_base);

    // ���������� �-�������
    // TODO: ������� ��� ��������� ��������
    // TODO: ����� ��������� boost
    double factorial(double k);

    // ��������� �������� ������� �� ������� k = 1, 2, 3 � ����� x ��� ��������
    // �������� t
    double FermiDirakFunction(double t, double x, double k);

    // ��������� �������� ������� �� ������� k = 1, 2, 3 � ����� x �� ����� 
    // ������� ��� x <= -0.1. N - ����� ������ � ����� ������� ��� ���������� 
    // �������� ��������
    double Gorner(double x, int N, double k);

    // ��������� �������� ������� �� ������� k = 1, 2, 3 � ����� x �� ������� 
    // ������-����������� � 5 ������ �� ����������� �����. 
    // (������� ������� �����) N = 5
    // TODO: std::function
    double FDGK5(double(*Ft)(double, double, double), double x, double T, double k, int N);

    // �������� �� ���������� ���������� ������ ������� FDGK5 
    double Richardson_mesh_refinement(double x, double t, double k, int N, int p);
} // namespace fdsfInteger