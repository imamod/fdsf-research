// Implementation of difference_scheme.h
//
#include "stdafx.h"
#include "matrix_helper.h"
#include "difference_scheme.h"

// ������� ��
double func_Fermi(double t, double x, double k) 
{
    // testfunction for positive X
    double u = pow(t,k)/(1 + exp(t-x));
    // ��� ������� ����������, ��� ����� �������
    //double u = pow(t, k) / (factorial(k)*(exp(x) + exp(t)));
    return u;
}

// ����� �������
double Gorner(double x, int N, double k)
{
    double exp_x = exp(x);
    double sum = 1.0 / pow(N, k + 1);

    for (int i = N - 1; i > 0; i--){
        sum = 1 / pow(i, k + 1) - exp_x*sum;
    }

    return sum;
}

// ���������� �-������� TODO: ������� ��� ��������� ��������
double factorial(double k)
{
    if (k < 0)
        return 0;
    if (k == 0)
        return 1;
    else
        return k*factorial(k - 1);
}

double regular_mesh_with_three_points(double(*Ft)(double, double, double), double x, double T, double k, int N)
{
    matrix_type::_vector gamma = { 5.0 / 18.0, 8.0 / 18.0 };
    matrix_type::_vector t(3);

    double U = 0;

    for (int n = 0; n < N; n++) {
        t.at(0) = T*(n + 0.5 - 0.5*sqrt(0.6)) / N;
        t.at(1) = T*(n + 0.5) / N;
        t.at(2) = T*(n + 0.5 + 0.5*sqrt(0.6)) / N;
        U = U + T*(Ft(t.at(1), x, k)*gamma.at(1) + gamma.at(0)*((Ft(t.at(0), x, k)) + Ft(t.at(2), x, k)));
    }

    U = U / N;

    return U;
}

// ����������� ����� � �������������� ������ ������-����������� (������� ������� �����) n = 5 ������������ �����
double regular_mesh_with_five_points(double(*Ft)(double, double, double), double x, double T, double k, int N)
{
    double gamma_1_5 = (322.0 - 13.0*sqrt(70)) / 1800.0;
    double gamma_2_4 = (322.0 + 13.0*sqrt(70)) / 1800.0;
    double gamma_3 = 64.0 / 225.0;
    matrix_type::_vector t(5);

    double U = 0;

    for (int n = N-1; n >=0 ; n--) {
        t.at(0) = T*(n + 0.5 - 0.5*sqrt((35 + 2 * sqrt(70)) / 63)) / N;
        t.at(1) = T*(n + 0.5 - 0.5*sqrt((35 - 2 * sqrt(70)) / 63)) / N;
        t.at(2) = T*(n + 0.5) / N;
        t.at(3) = T*(n + 0.5 + 0.5*sqrt((35 - 2 * sqrt(70)) / 63)) / N;
        t.at(4) = T*(n + 0.5 + 0.5*sqrt((35 + 2 * sqrt(70)) / 63)) / N;
        U = U + T*(Ft(t.at(2), x, k)*gamma_3 + gamma_1_5*((Ft(t.at(0), x, k)) + Ft(t.at(4), x, k)) + gamma_2_4*((Ft(t.at(1), x, k)) + Ft(t.at(3), x, k)));
    }

    U = U / N;

    return U;
}

// �������� �� ���������� �� �������-�������� ������
double Richardson_mesh_refinement(double x, double t, double k, int N, int p)
{
    double current_accuracy;
    double I_n, I_2n, I;

    //I_n = method_quasiregular_mesh_with_GK(&func_Fermi, x, k, N);
    //I_n = regular_mesh_with_three_points(&func_Fermi, x, t, k, N);
    I_n = regular_mesh_with_five_points(&func_Fermi, x, t, k, N);
    do {
        //I_2n = method_quasiregular_mesh_with_GK(&func_Fermi, x, k, N);
        //I_2n = regular_mesh_with_three_points(&func_Fermi, x, t, k, 2 * N);
        I_2n = regular_mesh_with_five_points(&func_Fermi, x, t, k, 2 * N);
        current_accuracy = (I_2n - I_n) / (pow(2.0, p) - 1);
        I = I_2n + current_accuracy;
        I_n = I_2n;
        N = 2 * N;
        //std::cout << "N = " << N << " R = " << current_accuracy << std::endl;
    } while (abs(current_accuracy) > epsilon); // ����������� �������� 10^-16

    //std::cout << std::endl << "    " << std::fixed << std::setprecision(std::numeric_limits<double>::digits10 + 2)
    //    << (I / log(2) - 1)*pow(10,10) << " " << log(2) << std::endl;
    return I;
}

// ���������������� ����� � �������������� ������ ������-����������� (3 ����)
double method_quasiregular_mesh_with_GK(double(*Ft)(double, double, double), double x, double k, int N)
{
    double m = 0.5, c = sqrt(3)*k*(1 + exp(-k));
    matrix_type::_vector ksi(3);
    matrix_type::_vector gamma = { 5.0 / 18.0, 8.0 / 18.0 };
    matrix_type::_vector t(3);
    matrix_type::_vector t_diff;

    double U = 0;

    for (int n = 0; n < N; n++) {
        ksi.at(0) = (n + 0.5 - 0.5*sqrt(0.6)) / N;
        ksi.at(1) = (n + 0.5) / N;
        ksi.at(2) = (n + 0.5 + 0.5*sqrt(0.6)) / N;
        for (int index = 0; index < 3; index++) {
            t.at(index) = c*ksi.at(index) / pow(1 - ksi.at(index)*ksi.at(index), m);
        }
        t_diff.push_back(c*(1 - (1 - 2 * m)*ksi.at(1)*ksi.at(1)) / pow(1 - ksi.at(1)*ksi.at(1), m + 1));
        U = U + (Ft(t.at(1), x, k)*gamma.at(1) + gamma.at(0)*((Ft(t.at(0), x, k)) + Ft(t.at(2), x, k)))*t_diff.at(n);
    }

    U = U / N;

    return U;
}

// ���������������� �����
double method_quasiregular_mesh(double(*Ft)(double, double, double), double x, double k, int N)
{
    double m = 0.5, c = sqrt(3)*k*(1 + exp(-k));
    matrix_type::_vector ksi;
    matrix_type::_vector t;
    matrix_type::_vector t_diff;

    double U = 0;

    for (int n = 0; n < N; n++) {
        ksi.push_back((n + 1 - 0.5) / N);
        t.push_back(c*ksi.at(n) / pow(1 - ksi.at(n)*ksi.at(n), m));
        t_diff.push_back(c*(1 - (1 - 2 * m)*ksi.at(n)*ksi.at(n)) / pow(1 - ksi.at(n)*ksi.at(n), m + 1));
        U = U + Ft(t.at(n), x, k)*t_diff.at(n);
    }

    U = U / N;

    return U;
}