#include "FDSFInteger.h"
#include <math.h>
#include <iomanip>
#include <limits>

mpfr::mpreal fdsf::get_T_max(mpfr::mpreal X, double k)
{
    mpfr::mpreal a1;
    if (k != 0) {
        a1 = mpfr::pow(fdsf::factorial(k + 1), -1 / k);
    }
    int i = 1;
    mpfr::mpreal y = mpfr::log(1 + mpfr::exp(X));
    mpfr::mpreal T = X - mpfr::log(epsilon), I_approximate;
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
            I_approximate = fdsf::factorial(k)*y*mpfr::pow(1 + a1*y, k);
        }
        // Итерационное вычисление Tmax
        T = X - mpfr::log(epsilon) + mpfr::log(mpfr::pow(T + 1, k) / I_approximate);
        i++;
    }

    return mpfr::sqrt(T);
}

// Новая схема Горнера для прецизионного вычисления y
mpfr::mpreal GornerSchemeForPrecesionY( mpfr::mpreal x, int N )
{
    //mpfr::mpreal alpha = x < epsilon ? exp(x) : exp(-x);
    const mpfr::mpreal alpha = mpfr::exp(x);
    const mpfr::mpreal z = alpha / (2 + alpha);
    mpfr::mpreal sum = 1.0 / (2 * N + 1);

    for (int i = N - 1; i >= 0; i--) {
        sum = 1.0 / (2 * i + 1.0) + z*z*sum;
    }

    //return x < epsilon ? 2*z*sum : x + 2*z*sum;
    return 2 * z*sum;
}

void fdsf::SetLinearTrigonometricGrid(std::vector<mpfr::mpreal> &y_base, 
                                      std::vector<mpfr::mpreal> &x_base,
                                      std::vector<mpfr::mpreal> &Y, 
                                      std::vector<mpfr::mpreal> &X, 
                                      const int N_base)
{
    int n_additional = 11;
    const mpfr::mpreal alpha = 2 / (2 + PI);
    const int base_nodes_count = 2 * N_base + 1;
    //Y = GornerSchemeForPrecesionY(N, X);
    // Задаются базовые узлы интерполяции
    for (int j = 1; j <= base_nodes_count; j++) {
        //y_base.push_back((1/2)*log(2)*(2 * alpha*j / (2 * N_base + 1) + (1 - alpha)*(1 - cos(PI*j / (2 * N_base + 1)))));
        y_base.push_back((mpfr::log(2)/2)*(2 * alpha*j / base_nodes_count +
                                          (1 - alpha)*(1 - mpfr::cos(PI*j / base_nodes_count))));
        x_base.push_back(mpfr::log(mpfr::exp(y_base.at(j - 1)) - 1));
    }

    Y.push_back(y_base.at(0) / n_additional);
    X.push_back(mpfr::log(mpfr::exp(Y.at(0)) - 1));

    for (int i = 1; i < n_additional; i++)
    {
        Y.push_back(Y.at(i - 1) + y_base.at(0) / n_additional);
        X.push_back(mpfr::log(mpfr::exp(Y.at(i)) - 1));
    }

    for (int index = 1; index < y_base.size(); index++) {
        for (int i = 0; i < n_additional; i++) {
            Y.push_back(Y.back() + (y_base.at(index) - y_base.at(index - 1)) / n_additional);
            X.push_back(mpfr::log(mpfr::exp(Y.back()) - 1));
        }
    }
}

mpfr::mpreal fdsf::FermiDirakFunction(mpfr::mpreal t, mpfr::mpreal x, mpfr::mpreal k)
{
//#if 0
    // new , зачин для полуцелых
    return 2 * mpfr::pow(t, 2 * k + 1) / (1 + mpfr::exp(t*t - x));
//#endif
#if 0
    return pow(t, k) / (factorial(k)*(exp(x) + exp(t)));
#endif
}

mpfr::mpreal fdsf::Gorner(mpfr::mpreal x, int N, mpfr::mpreal k)
{
    mpfr::mpreal exp_x = mpfr::exp(x);
    mpfr::mpreal sum = 1.0 / mpfr::pow(N, k + 1);

    for (int i = N - 1; i > 0; i--) {
        sum = 1 / mpfr::pow(i, k + 1) - exp_x*sum;
    }
#ifdef _DEBUG
    std::cout << "-----DEBUG-----" << std::endl;
    std::cout << sum << std::endl;
#endif
    return sum;
}

mpfr::mpreal fdsf::factorial(mpfr::mpreal k)
{
    if (k < 0)
        return 0;
    if (k == 0)
        return 1;
    else
        return k*factorial(k - 1);
}

mpfr::mpreal fdsf::FDGK5(mpfr::mpreal(*Ft)(mpfr::mpreal, mpfr::mpreal, mpfr::mpreal), 
                         const mpfr::mpreal x, 
                         const mpfr::mpreal T, 
                         const mpfr::mpreal k, 
                         const int N)
{
    // Веса формул Гаусса-Кристоффеля с N=5
    const mpfr::mpreal gamma_1_5 = (322.0 - 13.0*mpfr::sqrt(70)) / 1800.0;
    const mpfr::mpreal gamma_2_4 = (322.0 + 13.0*mpfr::sqrt(70)) / 1800.0;
    const mpfr::mpreal gamma_3 = 64.0 / 225.0;
    std::vector<mpfr::mpreal> t(5);

    const mpfr::mpreal half = 1.0/2.0;
    mpfr::mpreal U = 0;

    for (int n = N - 1; n >= 0; n--) {
        // Расчет дополнительных узлов
        t.at(0) = T*(n + half - (half)*mpfr::sqrt((35 + 2 * mpfr::sqrt(70)) / 63)) / N;
        t.at(1) = T*(n + half - (half)*mpfr::sqrt((35 - 2 * mpfr::sqrt(70)) / 63)) / N;
        t.at(2) = T*(n + half) / N;
        t.at(3) = T*(n + half + (half)*mpfr::sqrt((35 - 2 * mpfr::sqrt(70)) / 63)) / N;
        t.at(4) = T*(n + half + (half)*mpfr::sqrt((35 + 2 * mpfr::sqrt(70)) / 63)) / N;
        //
        U = U + T*(Ft(t.at(2), x, k)*gamma_3 + gamma_1_5*((Ft(t.at(0), x, k)) + Ft(t.at(4), x, k)) 
                                             + gamma_2_4*((Ft(t.at(1), x, k)) + Ft(t.at(3), x, k)));
    }

    U = U / N;

    return U;
}

// Сгущение по Ричардсону по сеточно-Гауссову методу
mpfr::mpreal fdsf::Richardson_mesh_refinement(mpfr::mpreal x, mpfr::mpreal t, mpfr::mpreal k)
{
    const int p = 10; // порядок аппроксимации формул ГК с 5ю узлами
    int N = 16; // начальное дробление сетки
    mpfr::mpreal current_accuracy, stop_criteria;
    mpfr::mpreal I_n, I_2n, I;
    // вспомогательное значение для расчета с повышенной разрядностью
    const mpfr::mpreal one = 1.0; 

    I_n = FDGK5(&FermiDirakFunction, x, t, k, N);
    do {
        I_2n = FDGK5(&FermiDirakFunction, x, t, k, 2 * N);
        current_accuracy = (I_2n - I_n) / (mpfr::pow(2.0, p) - 1);
        stop_criteria = (I_2n / I_n - one); // критерий останова подсчета
        I = I_2n;// +current_accuracy;
        I_n = I_2n;
        N = 2 * N;
        std::cout << "-----DEBUG-----" << std::endl;
        std::cout << I << std::endl;
    } while (mpfr::abs(stop_criteria) > 1e-20);
    //} while (abs(current_accuracy) > epsilon); // Фактическая точность 10^-16
#ifdef _DEBUG
    std::cout << I << std::endl;
#endif
    return I;
}

mpfr::mpreal fdsf::integer::FD_I1(mpfr::mpreal x)
{
    const mpfr::mpreal I_1_0 = fdsf::I_k_0[1];
    const mpfr::mpreal t = 60; 
    const int k = 1;

    mpfr::mpreal I_1_minus_x = fdsf::Richardson_mesh_refinement(-x, t, k);

    return x*x / 2 + 2*I_1_0 - I_1_minus_x;
}

mpfr::mpreal fdsf::integer::FD_I2(mpfr::mpreal x)
{
    const mpfr::mpreal I_1_0 = fdsf::I_k_0[1];
    const mpfr::mpreal t = 75; 
    const int k = 2;

    mpfr::mpreal I_2_minus_x = fdsf::Richardson_mesh_refinement(-x, t, k);

    return x*x*x / 3 + 4*x*I_1_0 + I_2_minus_x;
}

mpfr::mpreal fdsf::integer::FD_I3(mpfr::mpreal x)
{
    const mpfr::mpreal I_1_0 = fdsf::I_k_0[1], I_3_0 = fdsf::I_k_0[3];
    const mpfr::mpreal Tmax = 100; 
    const int k = 3;

    mpfr::mpreal I_3_minus_x = fdsf::Richardson_mesh_refinement(-x, Tmax, k);

    return x*x*x*x / 4 + 6*x*x*I_1_0 + 2*I_3_0 - I_3_minus_x;
}