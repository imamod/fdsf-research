#pragma once

#include <vector>
#include "mpreal.h"

namespace fdsf{

    // �������� �������� �������
    const int digits = 25;
    const mpfr::mpreal epsilon = pow(10,-digits);

    // �������� pi (����� ����� �� ���������� mpfr)
    //const mpfr::mpreal PI = 3.141592653589793238463;
    const mpfr::mpreal PI = mpfr::const_pi();

    // �������� Ik(0) ��� ������� k = 0,1,2,3
    const static mpfr::mpreal I_k_0[] = { log(2),
                                    PI*PI / 12.0, 
                                    1.8030853547393952,// �������� �� ������
                                    7.0*PI*PI*PI*PI / 120.0 };

    // ������ �������-������������������ ����� � ������� �����, ���� 10(?)
    // �������������� ����� ����� ������ ����� ������� �����.
    void SetLinearTrigonometricGrid(std::vector<mpfr::mpreal> &y_base, 
                                    std::vector<mpfr::mpreal> &x_base, 
                                    std::vector<mpfr::mpreal> &Y, 
                                    std::vector<mpfr::mpreal> &X, int N_base);

    // ���������� �-�������
    // TODO: ������� ��� ��������� ��������
    // TODO: ����� ��������� boost
    mpfr::mpreal factorial(mpfr::mpreal k);

    // ��������� �������� ������� �� ������� k = 1, 2, 3 � ����� x ��� ��������
    // �������� t. ������� ������������ � ���� t^k /(exp(x)+exp(t))
    mpfr::mpreal FermiDirakFunction(mpfr::mpreal t, mpfr::mpreal x, mpfr::mpreal k);

    // ��������� �������� ������� �� ������� k = 1, 2, 3 � ����� x �� ����� 
    // ������� ��� x <= -0.1. N - ����� ������ � ����� ������� ��� ���������� 
    // �������� ��������
    mpfr::mpreal Gorner(mpfr::mpreal x, int N, mpfr::mpreal k);

    // �������� ������ T_max ���������� x (� �� ������������� ��� ����)
    // NB: ��� ����� ����� ������ ��������������� �������
    mpfr::mpreal get_T_max(mpfr::mpreal X, double k);

    // ��������� �������� ������� �� ������� k = 1, 2, 3 � ����� x �� ������� 
    // ������-����������� � 5 ������ �� ����������� �����. 
    // �� ������, ��� ���������� ��������� �������� �� ������� -0.1<x<0 �����
    // ����� ��� k = 0 : Tmax = 50, 
    //           k = 1 : Tmax = 60,
    //           k = 2 : Tmax = 75, 
    //           k = 3 : Tmax = 100. 
    // TODO: std::function
    mpfr::mpreal FDGK5(mpfr::mpreal(*Ft)(mpfr::mpreal, mpfr::mpreal, mpfr::mpreal), 
                       mpfr::mpreal x, mpfr::mpreal T, mpfr::mpreal k, int N);

    // �������� �� ���������� ���������� ������ ������� FDGK5
    mpfr::mpreal Richardson_mesh_refinement(mpfr::mpreal x, mpfr::mpreal t, mpfr::mpreal k, int N);
    
    namespace integer {
        // ��������� �������� ������� �� ������� k=1 � ����� x
        mpfr::mpreal FD_I1(mpfr::mpreal x);

        // ��������� �������� ������� �� ������� k=2 � ����� x
        mpfr::mpreal FD_I2(mpfr::mpreal x);

        // ��������� �������� ������� �� ������� k=3 � ����� x
        mpfr::mpreal FD_I3(mpfr::mpreal x);
    } // namespace integer

} // namespace fdsf

mpfr::mpreal GornerSchemeForPrecesionY(mpfr::mpreal x, int N);