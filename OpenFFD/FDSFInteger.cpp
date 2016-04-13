#include "FDSFInteger.h"
#include <iomanip>
#include <limits>

using namespace fdsf;

cpp_dec_float_50 fdsf::get_T_max(cpp_dec_float_50 X, int k)
{
    cpp_dec_float_50 a1;
    if (k != 0) {
        a1 = pow(fdsf::factorial(k + 1), -1 / k);
    }
    int i = 1;
    cpp_dec_float_50 y = log(1 + exp(X));
    cpp_dec_float_50 T = X - log(epsilon), I_approximate;
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
static cpp_dec_float_50 GornerSchemeForPrecesionY(int N, cpp_dec_float_50 x)
{
    //cpp_dec_float_50 alpha = x < epsilon ? exp(x) : exp(-x);
    const cpp_dec_float_50 alpha = exp(x);
    const cpp_dec_float_50 z = alpha / (2 + alpha);
    cpp_dec_float_50 sum = 1.0 / (2 * N + 1);

    for (int i = N - 1; i >= 0; i--) {
        sum = 1.0 / (2 * i + 1.0) + z*z*sum;
    }

    //return x < epsilon ? 2*z*sum : x + 2*z*sum;
    return 2 * z*sum;
}

void fdsf::SetLinearTrigonometricGrid(std::vector<cpp_dec_float_50> &y_base, 
                                      std::vector<cpp_dec_float_50> &x_base,
                                      std::vector<cpp_dec_float_50> &Y, 
                                      std::vector<cpp_dec_float_50> &X, int N_base)
{
    int n_additional = 11;
    const cpp_dec_float_50 alpha = 2 / (2 + PI);
    const cpp_dec_float_50 one = cpp_dec_float_50(1);
    const cpp_dec_float_50 num2 = cpp_dec_float_50(2);
    cpp_dec_float_50 baseSize = cpp_dec_float_50(2 * N_base + 1);

    // Задаются базовые узлы интерполяции
    for (int j = 1; j <= baseSize; j++) {
        y_base.push_back(log(num2)/num2*(num2 * alpha*j / baseSize 
                                + (one - alpha)*(one - cos(PI*j / baseSize))));
        x_base.push_back(log(exp(y_base.at(j - 1)) - one));
    }

    Y.push_back(y_base.at(0) / n_additional);
    X.push_back(log(exp(Y.at(0)) - one));

    for (int i = 1; i < n_additional; i++)
    {
        Y.push_back(Y.at(i - 1) + y_base.at(0) / n_additional);
        X.push_back(log(exp(Y.at(i)) - one));
    }

    for (int index = 1; index < y_base.size(); index++) {
        for (int i = 0; i < n_additional; i++) {
            Y.push_back(Y.back() + (y_base.at(index) - y_base.at(index - 1)) / n_additional);
            X.push_back(log(exp(Y.back()) - one));
        }
    }
}

cpp_dec_float_50 fdsf::FermiDirakFunction(cpp_dec_float_50 t, cpp_dec_float_50 x, cpp_dec_float_50 k)
{
#if 0
    // new , зачин для полуцелых
    return 2 * pow(t, 2 * k + 1) / (1 + exp(t*t - x));
#endif
    return pow(t, k) / (factorial(k)*(exp(x) + exp(t)));
}

cpp_dec_float_50 fdsf::Gorner(cpp_dec_float_50 x, int N, cpp_dec_float_50 k)
{
    cpp_dec_float_50 exp_x = exp(x);
    cpp_dec_float_50 sum = 1.0 / pow(N, k + 1);

    for (int i = N - 1; i > 0; i--) {
        sum = 1 / pow(i, k + 1) - exp_x*sum;
    }

    return sum;
}

cpp_dec_float_50 fdsf::factorial(cpp_dec_float_50 k)
{
    if (k < 0)
        return 0;
    if (k == 0)
        return 1;
    else
        return k*factorial(k - 1);
}

cpp_dec_float_50 fdsf::FDGK5(cpp_dec_float_50(*Ft)(cpp_dec_float_50, cpp_dec_float_50, cpp_dec_float_50), 
                                                   cpp_dec_float_50 x, cpp_dec_float_50 T, cpp_dec_float_50 k, int N)
{
    //Определяем вспомогательные числа, чтобы не скатится до точности 16 знаков
    cpp_dec_float_50 half = cpp_dec_float_50(1) / 2;
    cpp_dec_float_50 num35 = cpp_dec_float_50(35);
    cpp_dec_float_50 num63 = cpp_dec_float_50(63);
    cpp_dec_float_50 num70 = cpp_dec_float_50(70);
    cpp_dec_float_50 num322 = cpp_dec_float_50(322);
    cpp_dec_float_50 num1800 = cpp_dec_float_50(1800);

    // Веса формул Гаусса-Кристоффеля с N=5
    const cpp_dec_float_50 gamma_1_5 = (num322 - cpp_dec_float_50(13)*sqrt(num70)) / num1800;
    const cpp_dec_float_50 gamma_2_4 = (num322 + cpp_dec_float_50(13)*sqrt(num70)) / num1800;
    const cpp_dec_float_50 gamma_3 = cpp_dec_float_50(64) / cpp_dec_float_50(225);
    std::vector<cpp_dec_float_50> t(5);

    cpp_dec_float_50 U = 0;

    for (int n = N - 1; n >= 0; n--) {
        // Расчет дополнительных узлов
        t.at(0) = T*(n + half - half*sqrt((num35 + 2 * sqrt(num70)) / num63)) / N;
        t.at(1) = T*(n + half - half*sqrt((num35 - 2 * sqrt(num70)) / num63)) / N;
        t.at(2) = T*(n + half) / N;
        t.at(3) = T*(n + half + half*sqrt((num35 - 2 * sqrt(num70)) / num63)) / N;
        t.at(4) = T*(n + half + half*sqrt((num35 + 2 * sqrt(num70)) / num63)) / N;
        //
        U = U + T*(Ft(t.at(2), x, k)*gamma_3 + gamma_1_5*((Ft(t.at(0), x, k)) + Ft(t.at(4), x, k)) 
                                             + gamma_2_4*((Ft(t.at(1), x, k)) + Ft(t.at(3), x, k)));
    }

    U = U / N;

    return U;
}

// Сгущение по Ричардсону по сеточно-Гауссову методу
cpp_dec_float_50 fdsf::Richardson_mesh_refinement(cpp_dec_float_50 x, cpp_dec_float_50 t, cpp_dec_float_50 k)
{
    int p = 10;
    int N = 16;
    cpp_dec_float_50 stop_criteria;
    cpp_dec_float_50 I_n, I_2n, I;

    I_n = FDGK5(&FermiDirakFunction, x, t, k, N);
    do {
        I_2n = FDGK5(&FermiDirakFunction, x, t, k, 2 * N);
        stop_criteria = (I_2n / I_n - 1);
        I = I_2n;
        I_n = I_2n;
        N = 2 * N;
    } while (abs(stop_criteria) > epsilon);
    
    return I;
}

cpp_dec_float_50 fdsf::integer::FD_I1(cpp_dec_float_50 x)
{
    const cpp_dec_float_50 I_1_0 = fdsf::I_k_0[1];
    const cpp_dec_float_50 t = 60; 
    const int k = 1;
  
    cpp_dec_float_50 I_1_minus_x = fdsf::Richardson_mesh_refinement(-x, t, k);

    return x*x / 2 + 2*I_1_0 - I_1_minus_x;
}

cpp_dec_float_50 fdsf::integer::FD_I2(cpp_dec_float_50 x)
{
    const cpp_dec_float_50 I_1_0 = fdsf::I_k_0[1];
    const cpp_dec_float_50 t = 75; 
    const int k = 2;

    cpp_dec_float_50 I_2_minus_x = fdsf::Richardson_mesh_refinement(-x, t, k);

    return x*x*x / 3 + 4*x*I_1_0 + I_2_minus_x;
}

cpp_dec_float_50 fdsf::integer::FD_I3(cpp_dec_float_50 x)
{
    const cpp_dec_float_50 I_1_0 = fdsf::I_k_0[1], I_3_0 = fdsf::I_k_0[3];
    const cpp_dec_float_50 Tmax = 100; 
    const int k = 3;

    cpp_dec_float_50 I_3_minus_x = fdsf::Richardson_mesh_refinement(-x, Tmax, k);

    return x*x*x*x / 4 + 6*x*x*I_1_0 + 2*I_3_0 - I_3_minus_x;
}