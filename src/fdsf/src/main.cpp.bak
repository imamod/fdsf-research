#include "Fdsf.h"
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
    BmpVector x = { -2, 2 };
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


// *****************************************************************************
// Функции работы с прецизионными аппроксимациями
// *****************************************************************************

void check_quadrature()
{
    const bmp_real k = bmp_real(1.0 / 2.0);
    bmp_real x = 30.0, I;
    bmp_real t = 0, a = 0;
    I = fdsf::richardson_method(x, t, k, a);
    std::cout << I << std::endl;
}


void check_z_value()
{
    BmpVector X, Y;
    const bmp_real k = bmp_real(1.0 / 2.0);
    bmp_real t = 0, a = 0, h = 0.1;
    BmpVector I;

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


BmpVector kuzmina_calculate_int(const bmp_real k) {
    BmpVector x0, X, y0, Y;
    BmpVector I_base, I_additional;
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
    BmpVector I_12_prec = kuzmina_calculate_int(bmp_real(1.0 / 2.0));
    BmpVector I_32_prec = kuzmina_calculate_int(bmp_real(3.0 / 2.0));
    BmpVector y;
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
    BmpVector alpha = { 58432.930, 21851.397, 1891.7921, 57.806497 };
    BmpVector betta = { 0.10730909, 0.0033951152 };
    BmpVector approximation;
    BmpVector delta;

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

// *****************************************************************************
// *****************************************************************************

int main()
{
    std::cout.precision(std::numeric_limits<bmp_real>::max_digits10);
#if 0
    printResultToFile(I_base, k, "I_base");
    printResultToFile(I_additional, k, "I_add");
    printResultToFile(y0, k, "y0");
    printResultToFile(Y, k, "Y");
#endif
    // *************************************************************************
    // Работа с прецизионными аппроксимациями полуцелых индексов
    // *************************************************************************
    //check_z_value();
    //check_quadrature();

    // Кузьмина
    //*******************************************
    //kuzminaCheck();
    //*******************************************
    // *************************************************************************
    // *************************************************************************
    //check();
    //control_point();
    //check_func_on_x();
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