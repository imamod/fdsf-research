#pragma once

#include "BasicTypes.h"

namespace fdsf {

    // �������� Ik(0) ��� ������� k = 0,1,2,3
    const static BmpReal I_k_0[] = { log(BmpReal(2)),
                                      PI*PI / 12.0, 
                                      1.8030853547393952,// �������� �� ������
                                      7.0*PI*PI*PI*PI / 120.0 
                                    };

    // ���������� �-�������
    BmpReal factorial(BmpReal k);

    // ��������� �������� ������� �� ������� k = 1, 2, 3 � ����� x ��� ��������
    // �������� t. ������� ������������ � ���� t^k /(exp(x)+exp(t))
    BmpReal fermi_dirak_integer(BmpReal t, BmpReal x, BmpReal k);
#if 0
    BmpReal fermi_dirak_half_integer(BmpReal ksi, BmpReal x, BmpReal k);
#endif

    // ��������� �������� ������� �� ������� k = 1, 2, 3 � ����� x �� ����� 
    // ������� ��� x <= -0.1. N - ����� ������ � ����� ������� ��� ���������� 
    // �������� ��������
    BmpReal Gorner(BmpReal x, BmpReal k);

    // �������� ������ T_max ���������� x (� �� ������������� ��� ����)
    // NB: ��� ����� ����� ������ ��������������� �������
    BmpReal get_T_max(BmpReal X, int k);

    // ��������� �������� ������� �� ������� k = 1, 2, 3 � ����� x �� ������� 
    // ������-����������� � 5 ������ �� ����������� �����. 
    // �� ������, ��� ���������� ��������� �������� �� ������� -0.1<x<0 �����
    // ����� ��� k = 0 : Tmax = 50, 
    //           k = 1 : Tmax = 60,
    //           k = 2 : Tmax = 75, 
    //           k = 3 : Tmax = 100. 
    // TODO: std::function
    BmpReal gauss_christoffel_method(BmpReal(*f)(BmpReal, BmpReal, BmpReal), 
                 BmpReal x, BmpReal T, BmpReal k, int N);

    // �������� �� ���������� ���������� ������ ������� gauss_christoffel_method
    BmpReal richardson_method(BmpReal x, BmpReal k, BmpReal t = 0);

    struct integration_segment_values {
        size_t n; // ������� ������� ��������������
        size_t N; // ����� ����� �������� ��������������
    };

    // for k = -3/2
    BmpReal fermi_dirak_m3half(BmpReal ksi, BmpReal x,
        BmpReal k, BmpReal a, integration_segment_values isv);

    // for others half-integer k
    BmpReal fermi_dirak_half_integer(BmpReal ksi, BmpReal x,
        BmpReal k, BmpReal a, integration_segment_values isv);

    // Formula Euler-Maclaurin
    BmpReal euler_maclaurin_method(BmpReal x, const BmpReal k, int N);

} // namespace fdsf
