#include "Fdsf.h"
#include <iomanip>
#include <limits>
#include <iostream>
#include <fstream>
#include <unordered_map>

namespace fdsf {

    bmp_real get_T_max(bmp_real X, int k)
    {
        bmp_real a1;
        if (k != 0) {
            a1 = pow(factorial(k + 1), -1 / k);
        }
        int i = 1;
        bmp_real y = log(1 + exp(X));
        bmp_real T = X - log(epsilon), I_approximate;
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
                I_approximate = factorial(k)*y*pow(1 + a1*y, k);
            }
            // Итерационное вычисление Tmax
            T = X - log(epsilon) + log(pow(T + 1, k) / I_approximate);
            i++;
        }

        return sqrt(T);
    }

    // Новая схема Горнера для прецизионного вычисления y
    static bmp_real GornerSchemeForPrecesionY(int N, bmp_real x)
    {
        //bmp_real alpha = x < epsilon ? exp(x) : exp(-x);
        const bmp_real alpha = exp(x);
        const bmp_real z = alpha / (2 + alpha);
        bmp_real sum = 1.0 / (2 * N + 1);

        for (int i = N - 1; i >= 0; i--) {
            sum = 1.0 / (2 * i + 1.0) + z*z*sum;
        }

        //return x < epsilon ? 2*z*sum : x + 2*z*sum;
        return 2 * z*sum;
    }

    void SetLinearTrigonometricGrid(std::vector<bmp_real> &y_base,
                                    std::vector<bmp_real> &x_base,
                                    std::vector<bmp_real> &Y,
                                    std::vector<bmp_real> &X, size_t N_base)
    {
        size_t n_additional = 11;
        const bmp_real alpha = 2 / (2 + PI);
        const bmp_real one = bmp_real(1);
        const bmp_real num2 = bmp_real(2); //if integer
        const bmp_real x_star = bmp_real(3);
        const bmp_real y_star = bmp_real(log(1 + exp(x_star))); // if half-integer
        bmp_real baseSize = bmp_real(2 * N_base + 1); // if integer || half-integer & !fixed a(N+1)
        //bmp_real baseSize = bmp_real(2 * N_base ); // if half-integer & fixed a(N+1)

        // Задаются базовые узлы интерполяции
        for (size_t j = 1; j <= baseSize; j++) {
            y_base.push_back(y_star / num2*(num2 * alpha*j / baseSize
                + (one - alpha)*(one - cos(PI*j / baseSize))));
            x_base.push_back(log(exp(y_base.at(j - 1)) - one));
        }

        Y.push_back(y_base.at(0) / n_additional);
        X.push_back(log(exp(Y.at(0)) - one));

        for (size_t i = 1; i < n_additional; i++)
        {
            Y.push_back(Y.at(i - 1) + y_base.at(0) / n_additional);
            X.push_back(log(exp(Y.at(i)) - one));
        }

        for (size_t index = 1; index < y_base.size(); index++) {
            for (size_t i = 0; i < n_additional; i++) {
                Y.push_back(Y.back() + (y_base.at(index) - y_base.at(index - 1)) / n_additional);
                X.push_back(log(exp(Y.back()) - one));
            }
        }
    }

    bmp_real fermi_dirak_integer(bmp_real t, 
                                 bmp_real x, 
                                 bmp_real k)
    {
        return pow(t, k) / (boost::math::tgamma(bmp_real(k)) * (exp(x) + exp(t)));
    }
#if 0
    bmp_real fermi_dirak_half_integer(bmp_real t, bmp_real x, bmp_real k)
    {
        // для полуцелых
        return 2 * pow(t, 2 * k + 1) / (1 + exp(t*t - x));
    }
#endif

    // Получение необходимого числа членов для схемы Горнера произвольной 
    // заданной точности
    static bmp_real get_N_for_Gorner(bmp_real x, bmp_real k)
    {
        return log(fdsf::epsilon) / x;
    }

    bmp_real Gorner(bmp_real x, bmp_real k)
    {
        //size_t N = ceil(get_N_for_Gorner(x, k));
        size_t N = (size_t)(get_N_for_Gorner(x, k)) + 1;
        bmp_real exp_x = exp(x);
        bmp_real sum = 1.0 / pow(N, k + 1);

        for (size_t i = N - 1; i > 0; i--) {
            sum = 1 / pow(i, k + 1) - exp_x * sum;
        }

        return sum * factorial(k) * exp(x);
    }

    // TODO: rename as gamma
    bmp_real factorial(bmp_real _k) {
        auto k = static_cast<double>(_k);
        std::unordered_map<double, bmp_real> SUPPORTED_HALFINTEGER_INDICES = {
            { -1.5, -2 * sqrt(PI) },
            { -0.5, sqrt(PI) },
            { 0.5, sqrt(PI) / 2 },
            { 1.5, 3 * sqrt(PI) / 4 },
            { 2.5, 15 * sqrt(PI) / 8 },
            { 3.5, 105 * sqrt(PI) / 16 }
        };

        auto it = SUPPORTED_HALFINTEGER_INDICES.find(k);
        if (it != SUPPORTED_HALFINTEGER_INDICES.end()) {
            return SUPPORTED_HALFINTEGER_INDICES[k];
        }
        if (!k) {
            return 1;
        }
        return k*factorial(k - 1);
    }

    bmp_real gauss_christoffel_method(bmp_real(*f)(bmp_real, bmp_real, bmp_real),
                                 bmp_real x, bmp_real T, bmp_real k, int N)
    {
        // Определяем вспомогательные числа, чтобы не скатится до точности 16 знаков
        bmp_real half = bmp_real(1.0) / 2;
        bmp_real num35 = bmp_real(35);
        bmp_real num63 = bmp_real(63);
        bmp_real num70 = bmp_real(70);
        bmp_real num322 = bmp_real(322);
        bmp_real num1800 = bmp_real(1800);

        // Веса формул Гаусса-Кристоффеля с N=5
        const bmp_real gamma_1_5 = (num322 - bmp_real(13)*sqrt(num70)) / num1800;
        const bmp_real gamma_2_4 = (num322 + bmp_real(13)*sqrt(num70)) / num1800;
        const bmp_real gamma_3 = bmp_real(64) / bmp_real(225);
        std::vector<bmp_real> t(5);

        bmp_real U = 0;

        for (int n = N - 1; n >= 0; n--) {
            // Расчет дополнительных узлов
            t[0] = T * (n + half - half * sqrt((num35 + 2 * sqrt(num70)) / num63)) / N;
            t[1] = T * (n + half - half * sqrt((num35 - 2 * sqrt(num70)) / num63)) / N;
            t[2] = T * (n + half) / N;
            t[3] = T * (n + half + half * sqrt((num35 - 2 * sqrt(num70)) / num63)) / N;
            t[4] = T * (n + half + half * sqrt((num35 + 2 * sqrt(num70)) / num63)) / N;
            //
            U = U + T * (f(t[2], x, k) * gamma_3 + gamma_1_5 * ((f(t[0], x, k)) + f(t[4], x, k))
                + gamma_2_4 * ((f(t[1], x, k)) + f(t[3], x, k)));
        }

        return U / N;
    }

    bool hasFractionalPart(bmp_real value) {
        return value - floor(value) > 0;
    }

    // Сгущение по Ричардсону по сеточно-Гауссову методу
    bmp_real richardson_method(bmp_real x,
                               bmp_real t,
                               bmp_real k,
                               bmp_real& a)
    {
        int N = 2;
        bmp_real stop_criteria;
        bmp_real I_n, I_2n;

        bool isHalfInteger = hasFractionalPart(k);
        if (isHalfInteger) {
            I_n = euler_maclaurin_method(x, k, N, a);
        }
        else {
            I_n = gauss_christoffel_method(&fermi_dirak_integer, x, t, k, N);
        }
        //std::ofstream fout;
        //fout.open("check_a_x.txt");
        //std::cout << "x = " << x << std::endl;
        //fout << (I_n) << std::fixed <<
        //    std::setprecision(std::numeric_limits<bmp_real>::max_digits10) << std::endl;
        //std::cout << "N = " << N << ": I = " << I_n << std::endl;
        do {
            if (isHalfInteger) {
                I_2n = euler_maclaurin_method(x, k, 2 * N, a);
                //I_2n += I_n / 2; // для трапеции если сетки по удвоенной
            }
            else {
                I_2n = gauss_christoffel_method(&fermi_dirak_integer, x, t, k, 2 * N);
            }
            //fout << I_2n << std::fixed <<
            //    std::setprecision(std::numeric_limits<bmp_real>::max_digits10) << std::endl;

            stop_criteria = (I_n / I_2n - 1); 
            //stop_criteria = (I_n / 0.52115038310799122 - 1); //k=-0.5
            I_n = I_2n;
            N = 2 * N;
            //std::cout << "N = " << N << ": d = "<< abs(stop_criteria) << std::endl; 
            //std::cout << "N = " << N << ": I = " << I_n << std::endl;
        } while (abs(stop_criteria) > 1e-11);
        //while (abs(stop_criteria) > epsilon*100);
        //std::cout << "N = " << N << ": I = ";
        //std::cout << I_2n << std::fixed <<
        //    std::setprecision(std::numeric_limits<bmp_real>::max_digits10) << std::endl;

        //fout.close();
        return I_n;
    }

}