#pragma once

#include <vector>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/detail/default_ops.hpp>
#include <boost/multiprecision/number.hpp>
#include <boost/math/special_functions/gamma.hpp>

namespace fdsf{
    using namespace boost::multiprecision;

    // �������� �������� �������
    const cpp_dec_float_50 epsilon = 1e-25;
    
    // �������� pi
    const cpp_dec_float_50 PI = boost::math::constants::pi<cpp_dec_float_50>();

    // �������� Ik(0) ��� ������� k = 0,1,2,3
    const static cpp_dec_float_50 I_k_0[] = { log(cpp_dec_float_50(2)),
                                              PI*PI / 12.0, 
                                              1.8030853547393952,// �������� �� ������
                                              7.0*PI*PI*PI*PI / 120.0 };

    // ������ �������-������������������ ����� � ������� �����, ���� 10(?)
    // �������������� ����� ����� ������ ����� ������� �����.
    void SetLinearTrigonometricGrid(std::vector<cpp_dec_float_50> &y_base, 
                                    std::vector<cpp_dec_float_50> &x_base, 
                                    std::vector<cpp_dec_float_50> &Y, 
                                    std::vector<cpp_dec_float_50> &X, int N_base);

    // ���������� �-�������
    // TODO: ������� ��� ��������� ��������
    // TODO: ����� ��������� boost
    cpp_dec_float_50 factorial(cpp_dec_float_50 k);

    // ��������� �������� ������� �� ������� k = 1, 2, 3 � ����� x ��� ��������
    // �������� t. ������� ������������ � ���� t^k /(exp(x)+exp(t))
    cpp_dec_float_50 FermiDirakFunction(cpp_dec_float_50 t, cpp_dec_float_50 x, cpp_dec_float_50 k);

    // ��������� �������� ������� �� ������� k = 1, 2, 3 � ����� x �� ����� 
    // ������� ��� x <= -0.1. N - ����� ������ � ����� ������� ��� ���������� 
    // �������� ��������
    cpp_dec_float_50 Gorner(cpp_dec_float_50 x, int N, cpp_dec_float_50 k);

    // �������� ������ T_max ���������� x (� �� ������������� ��� ����)
    // NB: ��� ����� ����� ������ ��������������� �������
    cpp_dec_float_50 get_T_max(cpp_dec_float_50 X, int k);

    // ��������� �������� ������� �� ������� k = 1, 2, 3 � ����� x �� ������� 
    // ������-����������� � 5 ������ �� ����������� �����. 
    // �� ������, ��� ���������� ��������� �������� �� ������� -0.1<x<0 �����
    // ����� ��� k = 0 : Tmax = 50, 
    //           k = 1 : Tmax = 60,
    //           k = 2 : Tmax = 75, 
    //           k = 3 : Tmax = 100. 
    // TODO: std::function
    cpp_dec_float_50 FDGK5(cpp_dec_float_50(*Ft)(cpp_dec_float_50, cpp_dec_float_50, cpp_dec_float_50), 
                 cpp_dec_float_50 x, cpp_dec_float_50 T, cpp_dec_float_50 k, int N);

    // �������� �� ���������� ���������� ������ ������� FDGK5
    cpp_dec_float_50 Richardson_mesh_refinement(cpp_dec_float_50 x, cpp_dec_float_50 t, cpp_dec_float_50 k);
    
    namespace integer {
        // ��������� �������� ������� �� ������� k=1 � ����� x
        cpp_dec_float_50 FD_I1(cpp_dec_float_50 x);

        // ��������� �������� ������� �� ������� k=2 � ����� x
        cpp_dec_float_50 FD_I2(cpp_dec_float_50 x);

        // ��������� �������� ������� �� ������� k=3 � ����� x
        cpp_dec_float_50 FD_I3(cpp_dec_float_50 x);
    } // namespace integer

} // namespace fdsf