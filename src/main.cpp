#include "FDSFInteger.h"
#include "FDSFHalf.h"
#include <boost/math/special_functions/gamma.hpp>

#include "MatrixUtils.h"
#include <iomanip>
#include <limits>
#include <fstream>
#include <sstream>
#include <iostream>

#include "ExponentialConvergence.h"
#include "../tefis/src/fdi/fdsf/fdsf.h"

#if 0
void test()
{
    using namespace boost::multiprecision;
    std::setprecision(std::numeric_limits<bmp_real>::max_digits10);
    // Operations at fixed precision and full numeric_limits support:
    cpp_dec_float_100 b = 2;
    std::cout << std::numeric_limits<cpp_dec_float_100>::digits << std::endl;
    // Note that digits10 is the same as digits, since we're base 10! :
    std::cout << std::numeric_limits<cpp_dec_float_100>::digits10 << std::endl;
    // We can use any C++ std lib function, lets print all the digits as well:
    std::cout << std::setprecision(std::numeric_limits<cpp_dec_float_100>::max_digits10)
    << log(b) << std::endl << log(bmp_real(2)) << std::endl; // print log(2)
                            // We can also use any function from Boost.Math:
    std::cout << boost::math::tgamma(b) << std::endl;
    // These even work when the argument is an expression template:
    std::cout << boost::math::tgamma(b * b) << std::endl;
    // And since we have an extended exponent range we can generate some really large 
    // numbers here (4.0238726007709377354370243e+2564):
    std::cout << boost::math::tgamma(cpp_dec_float_100(1000)) << std::endl;
    std::cout << boost::math::tgamma(cpp_dec_float_100(3.0/2.0)) << std::endl;

    bmp_real my_pi = boost::math::constants::pi<bmp_real>();
    std::cout << "my_pi = " << my_pi << std::endl;
    std::cout << "4atg1 = " << 4*boost::multiprecision::atan(bmp_real(1)) << std::endl;
    std::cout << "I1(0) = " << my_pi*my_pi/12 << std::endl;
}
#endif
#if 0
void check()
{
    std::cout.precision(std::numeric_limits<bmp_real>::max_digits10);
    bmp_real k = 1.0 / 2;
    bmp_real t = 0;
    bmp_real x0 = -1.0;
    bmp_real I_base = fdsf::richardson_method(x0, t, k);
    std::cout << "I_base: " << I_base << std::endl;
//#if 0
    bmp_real I_prec = fdsf::Gorner(x0, k);
    std::cout << "I_prec: " << I_prec << std::endl;
//#endif
}
#endif
// Сравнение схемы Горнера и метода трапеций
void compare()
{
    std::cout.precision(std::numeric_limits<bmp_real>::max_digits10);
    bmp_real k = 1.0 / 2;
    bmp_real x = -1;
    bmp_real a;
    bmp_real I_base = fdsf::richardson_method(x, 0, k, a);
    bmp_real I_prec = fdsf::Gorner(x, k);
    std::cout << I_base << std::endl;
    std::cout << I_prec << std::endl;
}

static bmp_real getN(bmp_real x, bmp_real k)
{
    bmp_real N0 = 1;
    bmp_real N = 1 + (log(epsilon) + (k+1)*log(N0))/x;

    while (true) {
        if (abs(N - N0) < 0.1) {
            break;
        }
        N0 = N;
        N = 1 + (log(epsilon) + (k + 1)*log(N0)) / x;
    }
    std::cout << N << std::endl;
    return N;
}

static bmp_real controlGorner(bmp_real x, bmp_real k)
{
    // size_t N = ceil(getN(x, k));
    size_t N = (size_t)(getN(x, k)) + 1;
    std::cout << N << std::endl;
    bmp_real exp_x = exp(x);
    bmp_real sum = 1.0 / pow(N, k + 1);

    for (size_t i = N - 1; i > 0; i--) {
        sum = 1 / pow(i, k + 1) - exp_x*sum;
    }

    return sum*factorial(k)*exp(x);
}

void control_point()
{
    std::cout.precision(std::numeric_limits<bmp_real>::max_digits10);
    bmp_real k = -1.0 / 2;
    bmp_real x = -1;

    std::cout << controlGorner( x, k) << std::endl; 
}


void check_func_on_x()
{
    std::cout.precision(std::numeric_limits<bmp_real>::max_digits10);
    bmp_real k = bmp_real(-3.0 / 2);
    bmp_real x = -10.0;
    std::ofstream fout;
    //fout.open("check_a_x.txt");
    fout.open("I_k_m15.txt");
    fout.precision(std::numeric_limits<bmp_real>::max_digits10);
    while (true)
    {
        if (x > 50.1) {
            break;
        }

        bmp_real a;
        bmp_real I_base = fdsf::richardson_method(x, 0, k, a);
        std::cout << "I_base: " << I_base << std::endl;
        //fout << "x = " << x << " a = " << a << std::endl;
        fout << I_base << std::endl;
        x += 0.1;
    }
    fout.close();
}

void probe_dots()
{
    std::cout.precision(std::numeric_limits<bmp_real>::max_digits10);
    bmp_real k = 7.0 / 2;
    bmp_real t = 0;
    std::vector<bmp_real> x = { -2, 2 };
    for (size_t i = 0; i < x.size(); i++) {
        bmp_real a;
        bmp_real I_base = fdsf::richardson_method(x[i], t, k, a);
        std::cout << I_base << std::endl;
    }
}

#if 0
void forPlot()
{
        int N = 2;
        bmp_real k = 1.0 / 2;
        bmp_real x = -1.0;
        bmp_real stop_criteria;
        bmp_real I_n, I_2n;

        I_n = euler_maclaurin_method(x, k, N);

        std::ofstream fout;
        fout.open("forPlot.txt");
        fout << (I_n) << std::fixed <<
            std::setprecision(std::numeric_limits<bmp_real>::max_digits10) << std::endl;
        std::cout << "N = " << N << ": I = " << I_n << std::endl;
        do {
            I_2n = euler_maclaurin_method(x, k, N++);
            fout << I_2n << std::fixed <<
                std::setprecision(std::numeric_limits<bmp_real>::max_digits10) << std::endl;

            stop_criteria = (I_2n / I_n - 1);
            I_n = I_2n;
            std::cout << "N = " << N << ": I = " << I_n << std::endl;
        } while (N < 65);
        //while (abs(stop_criteria) > epsilon * 100);

        fout.close();
}
#endif
// Вывод значений в файл
void printResultToFile(matrix_type::_vector x, bmp_real k, std::string varName)
{
    std::string fileName = { 0 };
    std::stringstream convertor;
    convertor << k;
    convertor >> fileName;
    fileName = varName + "_k" + fileName + ".txt";

    std::ofstream f_out;
    f_out.open(fileName);
    f_out.precision(std::numeric_limits<bmp_real>::max_digits10);
    for (size_t i = 0; i < x.size(); i++) {
        f_out << std::fixed << x[i] << std::endl;
    }

    f_out.close();
}

void GetValue_w(matrix_type::_vector &I_base, 
                matrix_type::_vector &I_additional,
                matrix_type::_vector y0, 
                matrix_type::_vector x0, 
                matrix_type::_vector Y, 
                matrix_type::_vector X, bmp_real k)
{
    bmp_real y;
    for (size_t i = 0; i < I_base.size(); i++) {
        y = log(1 + exp(x0[i]));
        I_base[i] = pow((I_base[i]*exp(x0[i]) / y0[i]), 1.0 / k);
    }

    for (size_t i = 0; i < I_additional.size(); i++) {
        y = log(1 + exp(X[i]));
        I_additional[i] = pow((I_additional[i]*exp(X[i]) / Y[i]), 1 / k);
    }
}
#if 0
// Для целых
static void computeIntegral(std::vector<bmp_real> x0, 
                            std::vector<bmp_real> X,
                            std::vector<bmp_real>& I_base, 
                            std::vector<bmp_real> &I_additional,
                            bmp_real x_div, bmp_real k)
{
    for (int i = 0; i < x0.size(); i++) {
        if (x0.at(i) > x_div) {
            bmp_real t;
            switch ((int)k) {
                case 1:
                    t = 60;
                    break;
                case 2:
                    t = 75;
                    break;
                case 3:
                    t = 100;
                    break;
                case 4:
                    t = 120;
                    break;
                default: break;
            }
            I_base.push_back(fdsf::richardson_method(x0.at(i), t, k));
        }
        else {
            I_base.push_back(fdsf::Gorner(x0.at(i), k));
            //I_base.push_back(GornerSchemeForPrecesionY( x0.at(i), N));
        }
        std::cout << "x0: " << x0.at(i) << " I_base: " << I_base.at(i) << std::endl;
    }

   // std::cout << "      I1(0) : " << fdsf::PI*fdsf::PI / 12 << std::endl;
    std::cout << "x_div = " << x_div<< std::endl;
    for (int i = 0; i < X.size(); i++) {
        if (X.at(i) > x_div) {
            //bmp_real t = fdsf::get_T_max(X.at(i), k);
            //bmp_real t = 60.0;
            //bmp_real t = 75.0;
            bmp_real t = 120;
            I_additional.push_back(fdsf::richardson_method(X.at(i), t, k));
        }
        else {
            I_additional.push_back(fdsf::Gorner(X.at(i), k));
        }
        std::cout << "X: " << X.at(i) << " I_add: " << I_additional.at(i) << std::endl;
        //<< std::fixed << std::setprecision(std::numeric_limits<bmp_real>::max_digits10)
    }

}
#endif
// *****************************************************************************
// Функции работы с прецизионными аппроксимациями
// *****************************************************************************
static void computeIntegral(std::vector<bmp_real> x0,
                            std::vector<bmp_real> X,
                            std::vector<bmp_real>& I_base,
                            std::vector<bmp_real> &I_additional,
                            bmp_real k)
{
    bmp_real t = 0, a = 0;
    for (size_t i = 0; i < x0.size(); i++) {
        I_base.push_back(fdsf::richardson_method(x0.at(i), t, k, a));
        //std::cout << "x0: " << x0.at(i) << " I_base: " << I_base.at(i) << std::endl;
    }

    for (size_t i = 0; i < X.size(); i++) {
        I_additional.push_back(fdsf::richardson_method(X.at(i), t, k, a));
        //std::cout << "X: " << X.at(i) << " I_add: " << I_additional.at(i) << std::endl;
        //<< std::fixed << std::setprecision(std::numeric_limits<bmp_real>::max_digits10)
    }

}

void SetLinearTrigonometricGridRight(std::vector<bmp_real> &y_base,
                                     std::vector<bmp_real> &x_base,
                                     std::vector<bmp_real> &Y,
                                     std::vector<bmp_real> &X, size_t N_base)
{
    size_t n_additional = 11;
    const bmp_real alpha = 2 / (2 + PI);
    const bmp_real one = bmp_real(1);
    const bmp_real num2 = bmp_real(2); //if integer
    const bmp_real x_star = bmp_real(10);
    const bmp_real y_star = bmp_real(log(1 + exp(x_star))); // if half-integer
    //bmp_real baseSize = bmp_real(2 * N_base + 1); // if integer || half-integer & !fixed a(N+1)
    //bmp_real baseSize = bmp_real(2 * N_base ); // if half-integer & fixed a(N+1)
    bmp_real baseSize = bmp_real( N_base); // if poly approximation

    const bmp_real y_star_inv = 1 / (y_star * y_star);
    //const bmp_real y_star_inv = 1 / y_star;
    //const bmp_real y_star_inv = 1 / pow(y_star, 0.5);
    //const bmp_real y_star_inv = 1 / pow(y_star, 0.25 );
    //const bmp_real y_star_inv = 1 / pow(y_star, 3.0 / 2);

    // Задаются базовые узлы интерполяции
    for (size_t j = 1; j <= baseSize; j++) {
        y_base.push_back(y_star_inv / num2*(num2 * alpha*j / baseSize
            + (one - alpha)*(one - cos(PI*j / baseSize))));
    }

    // Задаются дополнительные точки
    Y.push_back(y_base[0] / n_additional);

    for (size_t i = 1; i < n_additional; i++)
    {
        Y.push_back(Y[i - 1] + y_base[0] / n_additional);
    }

    for (size_t index = 1; index < y_base.size(); index++) {
        for (size_t i = 0; i < n_additional; i++) {
            Y.push_back(Y.back() + (y_base[index] - y_base[index - 1]) / n_additional);
        }
    }

    // Разворачиваем y
    std::reverse(y_base.begin(), y_base.end());
    for (size_t j = 0; j < baseSize; j++) {
        y_base[j] = 1.0 / pow(y_base[j], 0.5);
        //y_base[j] = 1.0 / y_base[j];
        //y_base[j] = 1.0 / (y_base[j] * y_base[j]);
        //y_base[j] = 1.0 / (pow(y_base[j], 4));
        //y_base[j] = 1.0 / (pow(y_base[j], 2.0 / 3));
        x_base.push_back(log(exp(y_base[j]) - one));
    }

    std::reverse(Y.begin(), Y.end());
    //std::cout << Y.size() << std::endl;
    for (size_t j = 0; j < Y.size(); j++) {
        Y[j] = 1.0 / pow(Y[j], 0.5);
        //Y[j] = 1.0 / Y[j]; 
        //Y[j] = 1.0 / (Y[j] * Y[j]);
        //Y[j] = 1.0 / pow(Y[j], 4);
        //Y[j] = 1.0 / pow(Y[j], 2.0 / 3);
        X.push_back(log(exp(Y[j]) - one));
    }
}


std::vector<bmp_real> calculate_series_part(const bmp_real& k,
                                            std::vector<bmp_real>& X)
{
    std::vector<bmp_real> coeff_A = { pow(PI, 2) / 6.0,
                                      pow(PI, 4) / 90.0,
                                      pow(PI, 6) / 945.0,
                                      pow(PI, 8) / 9450.0,
                                      pow(PI, 10) / 93555.0,
                                      691.0 * pow(PI, 12) / 638512875.0
                                     };

    std::vector<bmp_real> series_value;

    for (size_t i = 0; i < X.size(); i++) {
        bmp_real coeff_C = 1;
        bmp_real nom = k + 1;
        bmp_real series_sum = bmp_real(1.0);
        for (size_t j = 0; j < coeff_A.size(); j++) {
            coeff_C *= nom*(nom-1); // По асимптотической формуле парное добавление множителей, поэтому далее отнимаем 2
            series_sum += 2.0 * (1.0 - pow(2.0, 1.0 - 2*(j+1))) * pow(X[i], (-2.0)*(j + 1))*coeff_A[j]*coeff_C;
            //std::cout << "A(j) = " << coeff_A[j] << ": series_sum = " << series_sum << std::endl;
            std::cout << "C(" << j+1 << ") = " << 2.0 * (1.0 - pow(2.0, 1.0 - 2 * (j + 1))) * coeff_A[j] * coeff_C << std::endl;
            nom -= 2;
        }
        series_sum *= pow(X[i], k + 1) / (k + 1);
        series_value.push_back(series_sum);
    }

    std::cout << PI*PI*PI*PI*PI*PI / 945.0 << std::endl;

    return series_value;

}

static bmp_real get_assympt_value(bmp_real x, bmp_real k)
{
    std::vector<bmp_real> I_minus, I, series_part, X = {x};
    bmp_real t = 0, a = 0;
    series_part = calculate_series_part(k, X);
    I_minus.push_back(fdsf::richardson_method(-x, t, k, a));
    I.push_back(I_minus[0] + series_part[0]);
    //return I[0];
    return series_part[0];
}

static bmp_real get_series_value(bmp_real x, bmp_real k)
{
    bmp_real series_value = 0;
    auto N = log(fdsf::epsilon) / (x);

    for (size_t n = 1; n < N; ++n) {
        series_value += pow(-1.0, n - 1) * exp(n*x) / pow(n, k + 1);
        //std::cout << series_value << std::endl;
    }

    return factorial(k)*series_value;
}

void comp_kostya_and_precise()
{
    const bmp_real k = bmp_real(1.0 / 2.0);
    bmp_real x = bmp_real(-0.1), I, I_kostya, I_precise;
    bmp_real t = 0, a = 0;
#if 0
    I = fdsf::richardson_method(x, t, k, a);
    I_kostya = fdsf::fd_half(x);
    I_precise = get_series_value(x, k);
    std::cout << "x = -0.1" << std::endl;
    std::cout << "I_quadrature: " << I << std::endl;
    std::cout << "I_kostya: " << I_kostya << std::endl;
    std::cout << "I_precise: " << I_precise << std::endl;
    std::cout << "delta = " << I / I_precise - 1 << std::endl;
#endif

    x = 30.0;// +10E-8;
    std::cout << "x = " << x << std::endl;
    I = fdsf::richardson_method(x, t, k, a);
    //I_kostya = fdsf::fd_half(x);
    I_precise = get_assympt_value(x, k);
    std::cout << "I_quadrature: " << I << std::endl;
    //std::cout << "I_kostya: " << I_kostya << std::endl;
    std::cout << "I_precise: " << I_precise << std::endl;
    std::cout << "delta = " << I / I_precise - 1 << std::endl;
}

void check_negative_quadrature_values()
{
    const bmp_real k = bmp_real(1.0 / 2.0);
    bmp_real x = bmp_real(-0.1), I, I_precise;
    bmp_real t = 0, a = 0;
    I = fdsf::richardson_method(x, t, k, a);
    I_precise = get_series_value(x, k);
    std::cout << "I_quadrature: " << I << std::endl;
    std::cout << "I_precise: " << I_precise << std::endl;
}

void check_quadrature()
{
    const bmp_real k = bmp_real(1.0 / 2.0);
    bmp_real x = 30.0, I;
    bmp_real t = 0, a = 0;
    I = fdsf::richardson_method(x, t, k, a);
    std::cout << I << std::endl;
}

void calculate_asimpt_value()
{
    std::vector<bmp_real> X, Y;
    const bmp_real k = bmp_real(1.0 / 2.0);
    bmp_real t = 0, a = 0, h = 0.1;
    std::vector<bmp_real> I, I_minus, series_part;

    //проверка на идиота при х = 30
    X.push_back(30.0);
    //X.push_back(log(exp(Y[0]) - 1));
    series_part = calculate_series_part(k, X);
    I_minus.push_back(fdsf::richardson_method(-X[0], t, k, a));
    I.push_back(I_minus[0] + series_part[0]);
    I.push_back(series_part[0]);
#if 0
    Y.push_back(3.0);
    X.push_back(log(exp(Y[0]) - 1));
    size_t i = 1;
    while (true)
    {
        Y.push_back(Y[0] + i*h);
        X.push_back(log(exp(Y[i]) - 1));

        if (Y[i] > 30.0) {
            break;
        }

        i++;
    }

    series_part = calculate_series_part(k, X);
    for (size_t i = 0; i < X.size(); i++) {
        I_minus.push_back(fdsf::richardson_method(-X[i], t, k, a));
        I.push_back( I_minus[i] + series_part[i]);
    }
#endif
    printResultToFile(I, k, "Asimpt_check");
}

void calculate_k_half_integer()
{
    std::vector<bmp_real> x0, X, y0, Y;
    std::vector<bmp_real> I_base, I_additional;
    const bmp_real k = bmp_real(7.0 / 2.0);
    const size_t N_base = 5;
    // Расчет значения интеграла в базовых узлах
    fdsf::SetLinearTrigonometricGrid(y0, x0, Y, X, N_base);
    //SetLinearTrigonometricGridRight(y0, x0, Y, X, N_base);
    // Расчет интеграла
    computeIntegral(x0, X, I_base, I_additional, k);

    printResultToFile(I_base, k, "I_base");
    printResultToFile(I_additional, k, "I_add");
    printResultToFile(y0, k, "y0");
    printResultToFile(Y, k, "Y");
}

void check_z_value()
{
    std::vector<bmp_real> X, Y;
    const bmp_real k = bmp_real(1.0 / 2.0);
    bmp_real t = 0, a = 0, h = 0.1;
    std::vector<bmp_real> I;

    Y.push_back(3.0);
    X.push_back(log(exp(Y[0]) - 1));
    size_t i = 1;
    while (true)
    {
        Y.push_back(Y[0] + i*h);
        X.push_back(log(exp(Y[i]) - 1));

        if (Y[i] > 50.0) {
            break;
        }

        i++;
    }

    for (size_t i = 0; i < X.size(); i++) {
        I.push_back(fdsf::richardson_method(X[i], t, k, a));
    }

    printResultToFile(I, k, "z_check");
}


std::vector<bmp_real> kuzmina_calculate_int(const bmp_real k) {
    std::vector<bmp_real> x0, X, y0, Y;
    std::vector<bmp_real> I_base, I_additional;
    const int N_base = 5;
    // Расчет значения интеграла в базовых узлах
    //fdsf::SetLinearTrigonometricGrid(y0, x0, Y, X, N_base);
    SetLinearTrigonometricGridRight(y0, x0, Y, X, N_base);
    // Расчет интеграла
    computeIntegral(x0, X, I_base, I_additional, k);
    return I_additional;
}

void kuzminaCheck() {
    double k_index = 3.0 / 2;
    std::vector<bmp_real> I_12_prec = kuzmina_calculate_int(bmp_real(1.0 / 2.0));
    std::vector<bmp_real> I_32_prec = kuzmina_calculate_int(bmp_real(3.0 / 2.0));
    std::vector<bmp_real> y;
    std::cout << "got I_32" << std::endl;
    std::cout << I_12_prec.size() << std::endl;
    std::cout << I_32_prec.size() << std::endl;
    for (size_t i = 0; i < I_12_prec.size(); ++i) {
        y.push_back(sqrt(3.0 * I_12_prec[i] / 2));
        //std::cout << "y = " << y[i] << std::endl;
    }
    std::cout << y.size() << std::endl;
    //double y_star = 3.75;
    double alpha0 = 2.0 / 5;
    std::vector<bmp_real> alpha = { 58432.930, 21851.397, 1891.7921, 57.806497 };
    std::vector<bmp_real> betta = { 0.10730909, 0.0033951152 };
    std::vector<bmp_real> approximation;
    std::vector<bmp_real> delta;

    ///TODO: correct calculation of I32
    auto k = 2; // from formula
    for (size_t i = 0; i < I_12_prec.size(); i++) {
        std::cout << "in cycle" << std::endl;
        bmp_real num;
        bmp_real denom;
        for (size_t p = 0; p < alpha.size(); p++ ) {
            num = alpha[p] * pow(y[i], 8.0*p / 3);
            if (p < betta.size()) {
                denom = betta[p] * pow(y[i], 8.0*p / 3);
            }
        }
        auto drob = num / denom;
        approximation.push_back(alpha0 * pow(y[i], k) * pow(drob, 1.0/4));
        delta.push_back(approximation[i]/I_32_prec[i]);
    }

    for (auto i : delta) {
        std::cout << "approx for k = 3/2: " << i-1 << std::endl;
    }
}

void calculate_all_k() {
    std::vector<bmp_real> k = { -0.5, 0.5, 1.5, 2.5, 3.5 };
    bmp_real left_start = -5.0, right_start = 3;
    bmp_real left_end = 5.0, right_end = 30;
    bmp_real span = 0.1;
    std::vector<bmp_real> x_left, x_right;
    for (auto i = left_start; i < left_end; i += span) {
        x_left.push_back(i);
    }
    for (auto i = right_start; i < right_end; i += span) {
        x_right.push_back(i);
    }
    //std::cout << x_left.size() << std::endl;
    for (auto index : k) {
        std::vector<bmp_real> I_left;
        for (auto item : x_left) {
            bmp_real t=0, a;
            I_left.push_back(fdsf::richardson_method(item, t, index, a));
            printResultToFile(I_left, index, "I");
        }
        std::vector<bmp_real> I_right;
        for (auto item : x_right) {
            bmp_real t = 0, a;
            I_right.push_back(fdsf::richardson_method(item, t, index, a));
            printResultToFile(I_right, index, "I_right");
        }
        //std::cout << I.size() << std::endl;
    }

    std::vector<bmp_real> I_m32;
    for (auto item : x_left) {
        //I_m32.push_back(fdsf::fd_3mhalf(item));
    }
    printResultToFile(I_m32, -1.5, "I");

    std::vector<bmp_real> I_m32_right;
    for (auto item : x_right) {
        //I_m32_right.push_back(fdsf::fd_3mhalf(item));
    }
    printResultToFile(I_m32_right, -1.5, "I_right");
}

void chebyshevBaseNodes(std::vector<bmp_real> &y_base,
                        std::vector<bmp_real> &x_base,
                        std::vector<bmp_real> &Y,
                        std::vector<bmp_real> &X, int N_base) {
    const size_t n_additional = 11;
    const bmp_real one = bmp_real(1);
    const bmp_real num2 = bmp_real(2); //if integer
    const bmp_real x_star = bmp_real(3);
    const bmp_real y_star = bmp_real(log(1 + exp(x_star))); // if half-integer
    bmp_real baseSize = bmp_real(N_base); // if poly approximation

    const bmp_real y_star_inv = 1 / (y_star * y_star);

    // Задаются базовые узлы интерполяции
    for (size_t j = 1; j <= baseSize; j++) {
        y_base.push_back( y_star_inv*pow(sin(PI*j/(2*baseSize)), 2) );
        std::cout << y_base[j - 1] << " ";
    }

    // Задаются дополнительные точки
    Y.push_back(y_base[0] / n_additional);

    for (size_t i = 1; i < n_additional; i++) {
        Y.push_back(Y[i - 1] + y_base[0] / n_additional);
    }

    for (size_t index = 1; index < y_base.size(); index++) {
        for (size_t i = 0; i < n_additional; i++) {
            Y.push_back(Y.back() + (y_base[index] - y_base[index - 1]) / n_additional);
        }
    }

    // Разворачиваем y
    std::reverse(y_base.begin(), y_base.end());
    for (size_t j = 0; j < baseSize; j++) {
        y_base[j] = 1.0 / pow(y_base[j], 0.5);
        x_base.push_back(log(exp(y_base[j]) - one));
    }

    std::reverse(Y.begin(), Y.end());
    //std::cout << Y.size() << std::endl;
    for (size_t j = 0; j < Y.size(); j++) {
        Y[j] = 1.0 / pow(Y[j], 0.5);
        X.push_back(log(exp(Y[j]) - one));
    }
}


void test_chebyshev_base_nodes() {
    bmp_real k = 0.5;
    std::vector<bmp_real> x0, X, y0, Y;
    std::vector<bmp_real> I_base, I_additional;
    const size_t N_base = 1;
    // Расчет значения интеграла в базовых узлах
    chebyshevBaseNodes(y0, x0, Y, X, N_base);
    // Расчет интеграла
    computeIntegral(x0, X, I_base, I_additional, k);

    printResultToFile(I_base, k, "I_base");
    printResultToFile(I_additional, k, "I_add");
    printResultToFile(y0, k, "y0");
    printResultToFile(Y, k, "Y");

}

// *****************************************************************************
// *****************************************************************************

int main()
{
    std::cout.precision(std::numeric_limits<bmp_real>::max_digits10);
    std::vector<bmp_real> x0, X, y0, Y;
    std::vector<bmp_real> I_base, I_additional;
    bmp_real x_div = bmp_real(-0.1);
    //std::cout << "x_div = " << x_div << std::endl;
    //int N_gorner = 260, k = 1;
    //int N_gorner = 214, k = 2;
    //int N_gorner = 165, k = 3;
    const bmp_real k = bmp_real(1.0/2.0);
    const size_t N_base = 2;

#if 0

    // Расчет значения интеграла в базовых узлах
    fdsf::SetLinearTrigonometricGrid(y0, x0, Y, X, N_base);

    // Расчет схемы Горнера и подсчета интеграла на Гауссовой сетке
    computeIntegral(x0, X, I_base, I_additional, x_div, k);
    //GetValue_w(I_base, I_additional, y0, x0, Y, X, k);

    printResultToFile(I_base, k, "I_base");
    printResultToFile(I_additional, k, "I_add");
    printResultToFile(y0, k, "y0");
    printResultToFile(Y, k, "Y");
#endif
    // *************************************************************************
    // Работа с прецизионными аппроксимациями полуцелых индексов
    // *************************************************************************
    //calculate_k_half_integer();
    //check_z_value();
    //calculate_asimpt_value();
    //check_quadrature();
    //check_negative_quadrature_values();
    //comp_kostya_and_precise();

    //calculate_all_k();
    /**********************************
     * Проверка чебышевской сетки
     **********************************/
    test_chebyshev_base_nodes();

    //*****************************

    // Кузьмина
    //*******************************************
    //kuzminaCheck();
    //*******************************************
    // *************************************************************************
    // *************************************************************************
    //epc::checkTrapz(0, PI); // Для статьи о сверхстепенной сходимости
    //epc::checkTrapz(0, bmp_real(1.0));
    //check();
    //control_point();
    //check_func_on_x();
    //compare();
    //probe_dots();
    //forPlot();
#if 0
    std::cout << "Begining with matrixes" << std::endl;
    CMatrix object;

    matrix_type::_matrix A, A_inv;
    matrix_type::_vector B, x, delta_base(y0.size(), 0), delta_add(Y.size(), 0);
    matrix_type::_vector a, b;
    object.fill_matrix(N_base, I_base, y0, B, A);

    //A_inv = inverse(A);
    //object.find_coefficients(A_inv, B, a, b, N_base);
    //printResultToFile(a, k, "a"); printResultToFile(b, k, "b");

    //GetApproxomateValues(a, b, y0, Y, I_additional, I_base, delta_base, delta_add, N_base);
    //printResultToFile(delta_base, k, "delta_base"); printResultToFile(delta_add, k, "delta_add");
#endif

    getchar();
    return 0;
}