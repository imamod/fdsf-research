#include <iostream>
#include "FDSFInteger.h"
#include "src\matrix_helper.h"
#include <vector>
#include "mpreal.h"
#include <iomanip>
#include <limits>
#include <fstream>

// Вывод значений в файл
void printResultToFile(matrix_type::_vector x, mpfr::mpreal k, std::string varName)
{
    std::string fileName = { 0 };

    switch ((int)k)
    {
    case 1:
    case 2:
    case 3:
    case 4:
    {
        std::stringstream convertor;
        convertor << k;
        convertor >> fileName;
        fileName = varName + "_k" + fileName + ".txt";
        break;
    }
    default:
        break;
    }

    std::ofstream fout; // создаём объект класса ofstream для записи 
    fout.open(fileName);
    fout.precision(digits);
    for (int i = 0; i < x.size(); i++)
    {
        fout << std::fixed << x.at(i) << std::endl;
    }

    fout.close();
}

void GetValue_w(matrix_type::_vector &I_base, 
                matrix_type::_vector &I_additional,
                matrix_type::_vector y0, 
                matrix_type::_vector x0, 
                matrix_type::_vector Y, 
                matrix_type::_vector X, mpfr::mpreal k)
{
    mpfr::mpreal y;
    for (int i = 0; i < I_base.size(); i++) {
        y = mpfr::log(1 + mpfr::exp(x0.at(i)));
        I_base.at(i) = mpfr::pow((I_base.at(i)*mpfr::exp(x0.at(i)) / y0.at(i)), 1 / k);
    }

    for (int i = 0; i < I_additional.size(); i++) {
        y = mpfr::log(1 + mpfr::exp(X.at(i)));
        I_additional.at(i) = mpfr::pow((I_additional.at(i)*mpfr::exp(X.at(i)) / Y.at(i)), 1 / k);
    }
}

static void computeIntegral(std::vector<mpfr::mpreal> x0, 
                            std::vector<mpfr::mpreal> X,
                            std::vector<mpfr::mpreal>& I_base, 
                            std::vector<mpfr::mpreal> &I_additional,
                            mpfr::mpreal x_div, mpfr::mpreal N_gorner, mpfr::mpreal k)
{
    for (int i = 0; i < x0.size(); i++) {
        if (x0.at(i) > x_div) {
            mpfr::mpreal t = fdsf::get_T_max(X.at(i), k);
            //std::cout << "t = " << t << std::endl;
            //mpfr::mpreal t = 60.0;
            //mpfr::mpreal t = 75.0;
            //mpfr::mpreal t = 110.0;
            I_base.push_back(fdsf::Richardson_mesh_refinement(x0.at(i), t, k));
        }
        else {
            I_base.push_back(fdsf::Gorner(x0.at(i), N_gorner, k));
            //I_base.push_back(GornerSchemeForPrecesionY( x0.at(i), N));
        }
        std::cout << "x0: " << x0.at(i) << " I_base: " << I_base.at(i) << std::endl;
    }

    std::cout << "      I1(0) : " << mpfr::const_pi()*mpfr::const_pi() / 12 << std::endl;
    std::cout << "x_div = " << x_div<< std::endl;
    for (int i = 0; i < X.size(); i++) {
        if (X.at(i) > x_div) {
            mpfr::mpreal t = fdsf::get_T_max(X.at(i), k);
            //mpfr::mpreal t = 60.0;
            //mpfr::mpreal t = 75.0;
            //mpfr::mpreal t = 110.0;
            I_additional.push_back(fdsf::Richardson_mesh_refinement(X.at(i), t, k));
        }
        else {
            I_additional.push_back(fdsf::Gorner(X.at(i), N_gorner, k));
        }
        std::cout << "X: " << X.at(i) << " I_add: " << I_additional.at(i) << std::endl;
        //<< std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10)
    }

}

mpfr::mpreal get_N_for_Gorner(mpfr::mpreal X, double k)
{
    int i = 1;
    mpfr::mpreal Sold = 0, S;
  
    while (true)
    {
        // Итерационное вычисление S
        S = -mpfr::log(fdsf::epsilon*mpfr::pow(Sold + 1, k + 1))/abs(X);
        Sold = S;
        i++;
        if (abs(Sold - S) < 0.1) {
            break;
        }
    }
    std::cout << "Iteration count = " << i << std::endl;
    return S;
}

void testOnHilbertMatrix()
{
    const int N = 12;
    const mpfr::mpreal one = 1.0;
    matrix_type::_matrix H(N, matrix_type::_vector(N, 0));

    for (size_t n = 0; n < N; n++) {
        for (size_t m = 0; m < N; m++) {
            H[n][m] = one / ((m+1)+(n+1)-1);
        }
    }

    matrix_type::_matrix H_inv = inverse(H);
    std::ofstream fout;
    fout.precision(digits);
    fout.open("out\\HILBERTTEST.txt");
    for (size_t n = 0; n < N; n++) {
        for (size_t m = 0; m < N; m++) {
            fout << H_inv[n][m] << " ";
        }
        fout << "\n\r";
    }
    fout.close();
}

int main()
{
    // Required precision of computations in decimal digits
    // Play with it to check different precisions
    const int digits = 30;
    //const mpfr::mpreal epsilon = pow(10, -digits);
    // Setup default precision for all subsequent computations
    // MPFR accepts precision in bits - so we do the conversion 
    mpfr::mpreal::set_default_prec(mpfr::digits2bits(digits));
    const mpfr::mpreal pi = mpfr::const_pi();

    std::cout.precision(digits);
    std::cout << "pi = " << pi << std::endl;

    std::vector<mpfr::mpreal> x0, X, y0, Y;
    std::vector<mpfr::mpreal> I_base, I_additional;
    mpfr::mpreal x_div = -0.1;
    std::cout << "x_div = " << x_div << std::endl;
    //int N_gorner = 260, k = 1;
    //int N_gorner = 214, k = 2;
    //int N_gorner = 165, k = 1;
    const mpfr::mpreal N_gorner = 461; //S = get_N_for_Gorner(x_div, k);
    const mpfr::mpreal k = 1;
    const mpfr::mpreal N_base = 4;

    //mpfr::mpreal S = get_N_for_Gorner(x_div, k);
    //std::cout << "N_gorner = " << S << std::endl;
    //testOnHilbertMatrix();
//#if 0
    //std::cout << "      I0(0) : " << log(2) << std::endl;
    //std::cout << " mpfr I0(0) : " << mpfr::log(2) << std::endl;

    // Расчет значения интеграла в базовых узлах
    fdsf::SetLinearTrigonometricGrid(y0, x0, Y, X, N_base);

    // Расчет схемы Горнера и подсчета интеграла на Гауссовой сетке
    computeIntegral(x0, X, I_base, I_additional, x_div, N_gorner, k);
    std::cout << "      I1(0)_teor : " << mpfr::const_pi()*mpfr::const_pi() / 12 << std::endl;
    GetValue_w(I_base, I_additional, y0, x0, Y, X, k);

    printResultToFile(I_base, k, "I_base");
    printResultToFile(I_additional, k, "I_add");
    printResultToFile(y0, k, "y0");
    printResultToFile(Y, k, "Y");
//#endif
//#if 0
    std::cout << "Begining with matrixes" << std::endl;
    CMatrix object;

    matrix_type::_matrix A, A_inv;
    matrix_type::_vector B, x, delta_base(y0.size(), 0), delta_add(Y.size(), 0);
#if 0
    matrix_type::_vector a = { 0.3126028287472988,
                               0.0673008212829461,
                               0.0087798043423074,
                               0.0007222414330882,
                               0.0000295873218273 },
                               b = {
                               0.0626028287472659,
                               0.0238723363198067,
                               0.0010727527758408,
                               0.0000687107172921 };
#endif
    matrix_type::_vector a, b;
    object.fill_matrix(N_base, I_base, y0, B, A);

    A_inv = inverse(A);
    object.find_coefficients(A_inv, B, a, b, N_base);
    printResultToFile(a, k, "a"); printResultToFile(b, k, "b");
    GetApproxomateValues(a, b, y0, Y, I_additional, I_base, delta_base, delta_add, N_base);
    printResultToFile(delta_base, k, "delta_base"); printResultToFile(delta_add, k, "delta_add");
//#endif

    getchar();
    return 0;
}