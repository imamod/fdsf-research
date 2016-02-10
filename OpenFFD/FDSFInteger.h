#pragma once

#include <vector>

namespace fdsfInteger{

    // �������� �������� �������
    const double epsilon = 1e-17;

    // �������� pi (����� ����� �� ���������� mpfr)
    const double PI = 3.141592653589793238463;

    // ��������� �����
    void setLinTrigMesh(std::vector<double> &y0, std::vector<double> &x0, std::vector<double> &Y, std::vector<double> &X, int N_base);

    // ���������� �-�������
    // TODO: ������� ��� ��������� ��������
    // TODO: ����� ��������� boost
    double factorial(double k);

    // ������� ��
    double func_Fermi(double t, double x, double k);

    // ����� �������
    double Gorner(double x, int N, double k);

    // ����������� ����� � �������������� ������ ������-����������� (������� ������� �����) N = 5
    double regular_mesh_with_five_points(double(*Ft)(double, double, double), double x, double T, double k, int N);

    // �������� �� ���������� �� �������-�������� ������
    double Richardson_mesh_refinement(double x, double t, double k, int N, int p);
}; // namespace fdsfInteger