#include "FDSFInteger.h"
#include <math.h>
#include <iomanip>
#include <limits>

mpfr::mpreal fdsf::get_T_max(mpfr::mpreal X, mpfr::mpreal k)
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
    return T;
    //return mpfr::sqrt(T);
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
                                      const mpfr::mpreal N_base)
{
    int n_additional = 11;
    const mpfr::mpreal one = 1.0;
    const mpfr::mpreal alpha = 2*one / (2*one + PI);
    const mpfr::mpreal base_nodes_count = 2 * N_base + 1;
    const mpfr::mpreal number2 = 2;
    //Y = GornerSchemeForPrecesionY(N, X);
    // Задаются базовые узлы интерполяции
    for (mpfr::mpreal j = 1; j <= base_nodes_count; j++) {
        //y_base.push_back((1/2)*log(2)*(2 * alpha*j / (2 * N_base + 1) + (1 - alpha)*(1 - cos(PI*j / (2 * N_base + 1)))));
        y_base.push_back((mpfr::log(number2)/number2)*(number2 * alpha*j / base_nodes_count +
                                          (one - alpha)*(one - mpfr::cos(PI*j / base_nodes_count))));
        x_base.push_back(mpfr::log(mpfr::exp(y_base.at((int)j - 1)) - one));
    }

    Y.push_back(y_base.at(0) / n_additional);
    X.push_back(mpfr::log(mpfr::exp(Y.at(0)) - one));

    for (int i = 1; i < n_additional; i++)
    {
        Y.push_back(Y.at(i - 1) + y_base.at(0) / n_additional);
        X.push_back(mpfr::log(mpfr::exp(Y.at(i)) - one));
    }

    for (int index = 1; index < y_base.size(); index++) {
        for (int i = 0; i < n_additional; i++) {
            Y.push_back(Y.back() + (y_base.at(index) - y_base.at(index - 1)) / n_additional);
            X.push_back(mpfr::log(mpfr::exp(Y.back()) - one));
        }
    }
}

mpfr::mpreal fdsf::FermiDirakFunction(const mpfr::mpreal t, const mpfr::mpreal x, const mpfr::mpreal k)
{
#if 0
    // new , зачин для полуцелых
    return 2 * mpfr::pow(t, 2 * k + 1) / (1 + mpfr::exp(t*t - x));
#endif
//#if 0
    // представление для целых
    // TODO: разнести по namespace-ам
    return mpfr::pow(t, k) / (factorial(k)*(mpfr::exp(x) + mpfr::exp(t)));
//#endif
}

mpfr::mpreal fdsf::Gorner(const mpfr::mpreal x, const mpfr::mpreal N, const mpfr::mpreal k)
{
    const mpfr::mpreal exp_x = mpfr::exp(x);
    const mpfr::mpreal one = 1.0;
    mpfr::mpreal sum = one / mpfr::pow(N, k + 1);

    for (mpfr::mpreal i = N - 1; i > 0; i--) {
        sum = one / mpfr::pow(i, k + 1) - exp_x*sum;
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
                         const mpfr::mpreal N)
{
    // Веса формул Гаусса-Кристоффеля с N=5
    const mpfr::mpreal number70 = 70;
    const mpfr::mpreal number13 = 13;
    const mpfr::mpreal number322 = 322;
    const mpfr::mpreal number1800 = 1800;

    const mpfr::mpreal gamma_1_5 = (number322 - number13*mpfr::sqrt(number70)) / number1800;
    const mpfr::mpreal gamma_2_4 = (number322 + number13*mpfr::sqrt(number70)) / number1800;
    const mpfr::mpreal gamma_3 = 64.0 / 225.0;
    std::vector<mpfr::mpreal> t(5);

    const mpfr::mpreal half = 1.0/2.0;
    mpfr::mpreal U = 0;

    const mpfr::mpreal number2 = 2;
    const mpfr::mpreal number35 = 35;
    const mpfr::mpreal number63 = 63;

    for (mpfr::mpreal n = N - 1; n >= 0; n--) {
        // Расчет дополнительных узлов
        t.at(0) = T*(n + half - (half)*mpfr::sqrt((number35 + number2 * mpfr::sqrt(number70)) / number63)) / N;
        t.at(1) = T*(n + half - (half)*mpfr::sqrt((number35 - number2 * mpfr::sqrt(number70)) / number63)) / N;
        t.at(2) = T*(n + half) / N;
        t.at(3) = T*(n + half + (half)*mpfr::sqrt((number35 - number2 * mpfr::sqrt(number70)) / number63)) / N;
        t.at(4) = T*(n + half + (half)*mpfr::sqrt((number35 + number2 * mpfr::sqrt(number70)) / number63)) / N;
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
    mpfr::mpreal N = 16; // начальное дробление сетки
    mpfr::mpreal current_accuracy, stop_criteria;
    mpfr::mpreal I_n, I_2n, I;
#ifdef _DEBUG
    // вспомогательное значение для расчета с повышенной разрядностью
    const mpfr::mpreal one = 1.0; 
#endif
    I_n = FDGK5(&FermiDirakFunction, x, t, k, N);
    do {
        I_2n = FDGK5(&FermiDirakFunction, x, t, k, 2 * N);
        stop_criteria = (I_2n / I_n - 1); // критерий останова подсчета
        I = I_2n;
        I_n = I_2n;
        N = 2 * N;
#ifdef _DEBUG
        std::cout << "-----DEBUG-----" << std::endl;
        std::cout << I << std::endl;
#endif
    } while (mpfr::abs(stop_criteria) > fdsf::epsilon);
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