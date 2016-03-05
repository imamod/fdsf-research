#pragma once

#include <vector>

namespace fdsf{

    // �������� �������� �������
    const double epsilon = 1e-17;

    // �������� pi (����� ����� �� ���������� mpfr)
    const double PI = 3.141592653589793238463;

    // �������� Ik(0) ��� ������� k = 0,1,2,3
    const static double I_k_0[] = { log(2),
                                    PI*PI / 12.0, 
                                    1.8030853547393952,// �������� �� ������
                                    7.0*PI*PI*PI*PI / 120.0 };

    // ������ �������-������������������ ����� � ������� �����, ���� 10(?)
    // �������������� ����� ����� ������ ����� ������� �����.
    void SetLinearTrigonometricGrid(std::vector<double> &y_base, 
                                    std::vector<double> &x_base, 
                                    std::vector<double> &Y, 
                                    std::vector<double> &X, int N_base);

    // ���������� �-�������
    // TODO: ������� ��� ��������� ��������
    // TODO: ����� ��������� boost
    double factorial(double k);

    // ��������� �������� ������� �� ������� k = 1, 2, 3 � ����� x ��� ��������
    // �������� t. ������� ������������ � ���� t^k /(exp(x)+exp(t))
    double FermiDirakFunction(double t, double x, double k);

    // ��������� �������� ������� �� ������� k = 1, 2, 3 � ����� x �� ����� 
    // ������� ��� x <= -0.1. N - ����� ������ � ����� ������� ��� ���������� 
    // �������� ��������
    double Gorner(double x, int N, double k);

    // �������� ������ T_max ���������� x (� �� ������������� ��� ����)
    // NB: ��� ����� ����� ������ ��������������� �������
    double get_T_max(double X, int k);

    // ��������� �������� ������� �� ������� k = 1, 2, 3 � ����� x �� ������� 
    // ������-����������� � 5 ������ �� ����������� �����. 
    // �� ������, ��� ���������� ��������� �������� �� ������� -0.1<x<0 �����
    // ����� ��� k = 0 : Tmax = 50, 
    //           k = 1 : Tmax = 60,
    //           k = 2 : Tmax = 75, 
    //           k = 3 : Tmax = 100. 
    // TODO: std::function
    double FDGK5(double(*Ft)(double, double, double), 
                 double x, double T, double k, int N);

    // �������� �� ���������� ���������� ������ ������� FDGK5
    double Richardson_mesh_refinement(double x, double t, double k, int N);
    
    namespace integer {
        // ��������� �������� ������� �� ������� k=1 � ����� x
        double FD_I1(double x);

        // ��������� �������� ������� �� ������� k=2 � ����� x
        double FD_I2(double x);

        // ��������� �������� ������� �� ������� k=3 � ����� x
        double FD_I3(double x);
    } // namespace integer

} // namespace fdsf