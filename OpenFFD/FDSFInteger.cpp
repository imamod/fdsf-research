#include "FDSFInteger.h"
#include <math.h>
#include <iomanip>
#include <limits>

double fdsf::get_T_max(double X, int k)
{
    double a1;
    if (k != 0) {
        a1 = pow(fdsf::factorial(k + 1), -1 / k);
    }
    int i = 1;
    double y = log(1 + exp(X));
    double T = X - log(epsilon), I_approximate;
    while (true)
    {
        // 3-х итераций вполне достаточно для определения Tmax
        if (i > 3) {
            break;
        }
        if (k == 0) {
            I_approximate = y;
        }
        else {
            I_approximate = fdsf::factorial(k)*y*pow(1 + a1*y, k);
        }
        // Итерационное вычисление Tmax
        T = X - log(epsilon) + log(pow(T + 1, k) / I_approximate);
        i++;
    }

    return sqrt(T);
}

// Новая схема Горнера для прецизионного вычисления y
static double GornerSchemeForPrecesionY(double N, double x)
{
    //double alpha = x < epsilon ? exp(x) : exp(-x);
    const double alpha = exp(x);
    const double z = alpha / (2 + alpha);
    double sum = 1.0 / (2 * N + 1);

    for (int i = N - 1; i >= 0; i--) {
        sum = 1.0 / (2 * i + 1.0) + z*z*sum;
    }

    //return x < epsilon ? 2*z*sum : x + 2*z*sum;
    return 2 * z*sum;
}

void fdsf::SetLinearTrigonometricGrid(std::vector<double> &y_base, 
                                      std::vector<double> &x_base,
                                      std::vector<double> &Y, 
                                      std::vector<double> &X, int N_base)
{
    int n_additional = 11;
    const double alpha = 2 / (2 + PI);

    // Задаются базовые узлы интерполяции
    for (int j = 1; j <= 2 * N_base + 1; j++) {
        y_base.push_back(0.5*log(2)*(2 * alpha*j / (2 * N_base + 1) + (1 - alpha)*(1 - cos(PI*j / (2 * N_base + 1)))));
        x_base.push_back(log(exp(y_base.at(j - 1)) - 1));
    }

    Y.push_back(y_base.at(0) / n_additional);
    X.push_back(log(exp(Y.at(0)) - 1));

    for (int i = 1; i < n_additional; i++)
    {
        Y.push_back(Y.at(i - 1) + y_base.at(0) / n_additional);
        X.push_back(log(exp(Y.at(i)) - 1));
    }

    for (int index = 1; index < y_base.size(); index++) {
        for (int i = 0; i < n_additional; i++) {
            Y.push_back(Y.back() + (y_base.at(index) - y_base.at(index - 1)) / n_additional);
            X.push_back(log(exp(Y.back()) - 1));
        }
    }
}

double fdsf::FermiDirakFunction(double t, double x, double k)
{
#if 0
    // new , зачин для полуцелых
    double u = 2 * pow(t, 2 * k + 1) / (1 + exp(t*t - x));
#endif
    return pow(t, k) / (factorial(k)*(exp(x) + exp(t)));
}

double fdsf::Gorner(double x, int N, double k)
{
    double exp_x = exp(x);
    double sum = 1.0 / pow(N, k + 1);

    for (int i = N - 1; i > 0; i--) {
        sum = 1 / pow(i, k + 1) - exp_x*sum;
    }

    return sum;
}

double fdsf::factorial(double k)
{
    if (k < 0)
        return 0;
    if (k == 0)
        return 1;
    else
        return k*factorial(k - 1);
}

double fdsf::FDGK5(double(*Ft)(double, double, double), double x, double T, double k, int N)
{
    // Веса формул Гаусса-Кристоффеля с N=5
    const double gamma_1_5 = (322.0 - 13.0*sqrt(70)) / 1800.0;
    const double gamma_2_4 = (322.0 + 13.0*sqrt(70)) / 1800.0;
    const double gamma_3 = 64.0 / 225.0;
    std::vector<double> t(5);

    double U = 0;

    for (int n = N - 1; n >= 0; n--) {
        // Расчет дополнительных узлов
        t.at(0) = T*(n + 0.5 - 0.5*sqrt((35 + 2 * sqrt(70)) / 63)) / N;
        t.at(1) = T*(n + 0.5 - 0.5*sqrt((35 - 2 * sqrt(70)) / 63)) / N;
        t.at(2) = T*(n + 0.5) / N;
        t.at(3) = T*(n + 0.5 + 0.5*sqrt((35 - 2 * sqrt(70)) / 63)) / N;
        t.at(4) = T*(n + 0.5 + 0.5*sqrt((35 + 2 * sqrt(70)) / 63)) / N;
        //
        U = U + T*(Ft(t.at(2), x, k)*gamma_3 + gamma_1_5*((Ft(t.at(0), x, k)) + Ft(t.at(4), x, k)) 
                                             + gamma_2_4*((Ft(t.at(1), x, k)) + Ft(t.at(3), x, k)));
    }

    U = U / N;

    return U;
}

// Сгущение по Ричардсону по сеточно-Гауссову методу
// TODO: пересмотреть критерий остановки
double fdsf::Richardson_mesh_refinement(double x, double t, double k)
{
    int p = 10;
    int N = 16;
    double current_accuracy;
    double I_n, I_2n, I;

    I_n = FDGK5(&FermiDirakFunction, x, t, k, N);
    do {
        I_2n = FDGK5(&FermiDirakFunction, x, t, k, 2 * N);
        current_accuracy = (I_2n - I_n) / (pow(2.0, p) - 1);
        I = I_2n + current_accuracy;
        I_n = I_2n;
        N = 2 * N;
    } while (abs(current_accuracy) > epsilon); // Фактическая точность 10^-16

    return I;
}

double fdsf::integer::FD_I1(double x)
{
    const double I_1_0 = fdsf::I_k_0[1];
    const double t = 60; 
    const int k = 1;
  
    double I_1_minus_x = fdsf::Richardson_mesh_refinement(-x, t, k);

    return x*x / 2 + 2*I_1_0 - I_1_minus_x;
}

double fdsf::integer::FD_I2(double x)
{
    const double I_1_0 = fdsf::I_k_0[1];
    const double t = 75; 
    const int k = 2;

    double I_2_minus_x = fdsf::Richardson_mesh_refinement(-x, t, k);

    return x*x*x / 3 + 4*x*I_1_0 + I_2_minus_x;
}

double fdsf::integer::FD_I3(double x)
{
    const double I_1_0 = fdsf::I_k_0[1], I_3_0 = fdsf::I_k_0[3];
    const double Tmax = 100; 
    const int k = 3;

    double I_3_minus_x = fdsf::Richardson_mesh_refinement(-x, Tmax, k);

    return x*x*x*x / 4 + 6*x*x*I_1_0 + 2*I_3_0 - I_3_minus_x;
}