#include "Fdsf-legacy.h"
#include "Constants.h"
#include "Gamma.h"
#include <iomanip>
#include <limits>
#include <iostream>
#include <fstream>

namespace fdsf {

    BmpReal get_T_max(BmpReal X, int k) {
        BmpReal a1;
        if (k != 0) {
            a1 = pow(factorial(k + 1), -1 / k);
        }
        int i = 1;
        BmpReal y = log(1 + exp(X));
        BmpReal T = X - log(epsilon), I_approximate;
        while (true) {
            // 3-х итераций вполне достаточно для определения Tmax
            if (i > 3) {
                break;
            }

            I_approximate = !k ? y : factorial(k)*y*pow(1 + a1*y, k);

            // Итерационное вычисление Tmax
            T = X - log(epsilon) + log(pow(T + 1, k) / I_approximate);
            i++;
        }

        return sqrt(T);
    }

    // Новая схема Горнера для прецизионного вычисления y
    static BmpReal GornerSchemeForPrecesionY(int N, BmpReal x) {
        //BmpReal alpha = x < epsilon ? exp(x) : exp(-x);
        const BmpReal alpha = exp(x);
        const BmpReal z = alpha / (2 + alpha);
        BmpReal sum = 1.0 / (2 * N + 1);

        for (int i = N - 1; i >= 0; i--) {
            sum = 1.0 / (2 * i + 1.0) + z*z*sum;
        }

        //return x < epsilon ? 2*z*sum : x + 2*z*sum;
        return 2 * z*sum;
    }

    BmpReal fermi_dirak_integer(BmpReal t, BmpReal x, BmpReal k) {
        return pow(t, k) / (boost::math::tgamma(BmpReal(k)) * (exp(x) + exp(t)));
    }

    // Получение необходимого числа членов для схемы Горнера произвольной 
    // заданной точности
    static BmpReal get_N_for_Gorner(BmpReal x, BmpReal k) {
        return log(fdsf::epsilon) / x;
    }

    BmpReal Gorner(BmpReal x, BmpReal k) {
        //size_t N = ceil(get_N_for_Gorner(x, k));
        size_t N = (size_t)(get_N_for_Gorner(x, k)) + 1;
        BmpReal exp_x = exp(x);
        BmpReal sum = 1.0 / pow(N, k + 1);

        for (size_t i = N - 1; i > 0; i--) {
            sum = 1 / pow(i, k + 1) - exp_x * sum;
        }

        return sum * factorial(k) * exp(x);
    }

    BmpReal gauss_christoffel_method(BmpReal(*f)(BmpReal, BmpReal, BmpReal),
                                 BmpReal x, BmpReal T, BmpReal k, int N) {
        // Определяем вспомогательные числа, чтобы не скатится до точности 16 знаков
        const BmpReal half = BmpReal(1.0) / 2;
        const BmpReal num35 = BmpReal(35);
        const BmpReal num63 = BmpReal(63);
        const BmpReal num70 = BmpReal(70);
        const BmpReal num322 = BmpReal(322);
        const BmpReal num1800 = BmpReal(1800);

        // Веса формул Гаусса-Кристоффеля с N=5
        const BmpReal gamma_1_5 = (num322 - BmpReal(13)*sqrt(num70)) / num1800;
        const BmpReal gamma_2_4 = (num322 + BmpReal(13)*sqrt(num70)) / num1800;
        const BmpReal gamma_3 = BmpReal(64) / BmpReal(225);
        BmpVector t(5);

        BmpReal U = 0;

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

    bool hasFractionalPart(BmpReal value) {
        return value - floor(value) > 0;
    }

    // Сгущение по Ричардсону по сеточно-Гауссову методу
    BmpReal richardson_method(BmpReal x, BmpReal k, BmpReal t) {
        int N = 2;
        BmpReal stop_criteria;
        BmpReal I_n, I_2n;

        bool isHalfInteger = hasFractionalPart(k);
        I_n = isHalfInteger ? euler_maclaurin_method(x, k, N) :
                              gauss_christoffel_method(&fermi_dirak_integer, x, t, k, N);
        do {
            if (isHalfInteger) {
                I_2n = euler_maclaurin_method(x, k, 2 * N);
                //I_2n += I_n / 2; // для трапеции если сетки по удвоенной
            } else {
                I_2n = gauss_christoffel_method(&fermi_dirak_integer, x, t, k, 2 * N);
            }

            stop_criteria = (I_n / I_2n - 1); 
            I_n = I_2n;
            N = 2 * N;
            if (N >= 2048) {
                std::cout << "N = " << N << ": I = " << I_2n;
            }
        } while (abs(stop_criteria) > 1e-11);
        //while (abs(stop_criteria) > epsilon*100);
        return I_n;
    }

}
