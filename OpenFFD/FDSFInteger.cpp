#include "FDSFInteger.h"
#include <math.h>
#include <iomanip>
#include <limits>

void fdsfInteger::setLinTrigMesh(std::vector<double> &y0, std::vector<double> &x0, std::vector<double> &Y, std::vector<double> &X, int N_base)
{
    int n_additional = 11;
    double alpha = 2 / (2 + PI);

    // Задается линейно-тригонометрическая сетка
    for (int j = 1; j <= 2 * N_base + 1; j++) {
        y0.push_back(0.5*log(2)*(2 * alpha*j / (2 * N_base + 1) + (1 - alpha)*(1 - cos(PI*j / (2 * N_base + 1)))));
        x0.push_back(log(exp(y0.at(j - 1)) - 1));
    }

    Y.push_back(y0.at(0) / n_additional);
    X.push_back(log(exp(Y.at(0)) - 1));

    for (int i = 1; i < n_additional; i++)
    {
        Y.push_back(Y.at(i - 1) + y0.at(0) / n_additional);
        X.push_back(log(exp(Y.at(i)) - 1));
    }

    for (int index = 1; index < y0.size(); index++) {
        for (int i = 0; i < n_additional; i++) {
            Y.push_back(Y.back() + (y0.at(index) - y0.at(index - 1)) / n_additional);
            X.push_back(log(exp(Y.back()) - 1));
        }
    }

}
// Функция ФД
double fdsfInteger::func_Fermi(double t, double x, double k)
{
    // testfunction for positive X
    //double u = pow(t,k)/(1 + exp(t-x));
    // Для расчета интегралом, Для схемы Горнера
    double u = pow(t, k) / (factorial(k)*(exp(x) + exp(t)));
    return u;
}

// Схема Горнера
double fdsfInteger::Gorner(double x, int N, double k)
{
    double exp_x = exp(x);
    double sum = 1.0 / pow(N, k + 1);

    for (int i = N - 1; i > 0; i--){
        sum = 1 / pow(i, k + 1) - exp_x*sum;
    }

    return sum;
}

// Вычисление Г-функции TODO: сделать для полуцелых индексов
double fdsfInteger::factorial(double k)
{
    if (k < 0)
        return 0;
    if (k == 0)
        return 1;
    else
        return k*factorial(k - 1);
}

// Равномерная сетка с использованием формул Гаусса-Кристофелля (сеточно Гауссов метод) n = 5 пятиточечная схема
double fdsfInteger::regular_mesh_with_five_points(double(*Ft)(double, double, double), double x, double T, double k, int N)
{
    double gamma_1_5 = (322.0 - 13.0*sqrt(70)) / 1800.0;
    double gamma_2_4 = (322.0 + 13.0*sqrt(70)) / 1800.0;
    double gamma_3 = 64.0 / 225.0;
    std::vector<double> t(5);

    double U = 0;

    for (int n = N - 1; n >= 0; n--) {
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

// Сгущение по Ричардсону по сеточно-Гауссову методу
double fdsfInteger::Richardson_mesh_refinement(double x, double t, double k, int N, int p)
{
    double current_accuracy;
    double I_n, I_2n, I;

    I_n = regular_mesh_with_five_points(&func_Fermi, x, t, k, N);
    do {
        I_2n = regular_mesh_with_five_points(&func_Fermi, x, t, k, 2 * N);
        current_accuracy = (I_2n - I_n) / (pow(2.0, p) - 1);
        I = I_2n + current_accuracy;
        I_n = I_2n;
        N = 2 * N;

    } while (abs(current_accuracy) > epsilon); // Фактическая точность 10^-16

    return I;
}