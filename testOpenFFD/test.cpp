#include "FDSFInteger.h"
#include <boost/math/special_functions/gamma.hpp>

#include "matrix_helper.h"
#include <iomanip>
#include <limits>
#include <fstream>
#include <sstream>
#include <iostream>

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

void check()
{
    std::cout.precision(std::numeric_limits<bmp_real>::max_digits10);
    bmp_real k = 1.0 / 2;
    bmp_real t = 0;
    bmp_real x0 = -1.0;
    bmp_real I_base = fdsf::Richardson_mesh_refinement(x0, t, k);
    std::cout << "I_base: " << I_base << std::endl;
//#if 0
    bmp_real I_prec = fdsf::Gorner(x0, k);
    std::cout << "I_prec: " << I_prec << std::endl;
//#endif
}

// Вывод значений в файл
void printResultToFile(matrix_type::_vector x, bmp_real k, std::string varName)
{
    std::string fileName = { 0 };
    std::stringstream convertor;
    convertor << k;
    convertor >> fileName;
    fileName = varName + "_k" + fileName + ".txt";

    std::ofstream fout; // создаём объект класса ofstream для записи 
    fout.open(fileName);
    fout.precision(std::numeric_limits<bmp_real>::max_digits10);
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
                matrix_type::_vector X, bmp_real k)
{
    bmp_real y;
    for (int i = 0; i < I_base.size(); i++) {
        y = log(1 + exp(x0.at(i)));
        I_base.at(i) = pow((I_base.at(i)*exp(x0.at(i)) / y0.at(i)), 1 / k);
    }

    for (int i = 0; i < I_additional.size(); i++) {
        y = log(1 + exp(X.at(i)));
        I_additional.at(i) = pow((I_additional.at(i)*exp(X.at(i)) / Y.at(i)), 1 / k);
    }
}

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
            I_base.push_back(fdsf::Richardson_mesh_refinement(x0.at(i), t, k));
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
            I_additional.push_back(fdsf::Richardson_mesh_refinement(X.at(i), t, k));
        }
        else {
            I_additional.push_back(fdsf::Gorner(X.at(i), k));
        }
        std::cout << "X: " << X.at(i) << " I_add: " << I_additional.at(i) << std::endl;
        //<< std::fixed << std::setprecision(std::numeric_limits<bmp_real>::max_digits10)
    }

}

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
    const int N_base = 4;

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

    check();

#if 0
    std::cout << "Begining with matrixes" << std::endl;
    CMatrix object;

    matrix_type::_matrix A, A_inv;
    matrix_type::_vector B, x, delta_base(y0.size(), 0), delta_add(Y.size(), 0);
    matrix_type::_vector a, b;
    //object.fill_matrix(N_base, I_base, y0, B, A);

    //A_inv = inverse(A);
    //object.find_coefficients(A_inv, B, a, b, N_base);
    //printResultToFile(a, k, "a"); printResultToFile(b, k, "b");

    //GetApproxomateValues(a, b, y0, Y, I_additional, I_base, delta_base, delta_add, N_base);
    //printResultToFile(delta_base, k, "delta_base"); printResultToFile(delta_add, k, "delta_add");
#endif

    getchar();
    return 0;
}