#include "FDSFInteger.h"
#include <boost/math/special_functions/gamma.hpp>

#include "src\matrix_helper.h"
#include <iomanip>
#include <limits>
#include <fstream>
#include <sstream>
#include <iostream>

void test()
{
    using namespace boost::multiprecision;
    std::setprecision(std::numeric_limits<cpp_dec_float_50>::max_digits10);
    // Operations at fixed precision and full numeric_limits support:
    cpp_dec_float_100 b = 2;
    std::cout << std::numeric_limits<cpp_dec_float_100>::digits << std::endl;
    // Note that digits10 is the same as digits, since we're base 10! :
    std::cout << std::numeric_limits<cpp_dec_float_100>::digits10 << std::endl;
    // We can use any C++ std lib function, lets print all the digits as well:
    std::cout << std::setprecision(std::numeric_limits<cpp_dec_float_100>::max_digits10)
    << log(b) << std::endl << log(cpp_dec_float_50(2)) << std::endl; // print log(2)
                            // We can also use any function from Boost.Math:
    std::cout << boost::math::tgamma(b) << std::endl;
    // These even work when the argument is an expression template:
    std::cout << boost::math::tgamma(b * b) << std::endl;
    // And since we have an extended exponent range we can generate some really large 
    // numbers here (4.0238726007709377354370243e+2564):
    std::cout << boost::math::tgamma(cpp_dec_float_100(1000)) << std::endl;
    std::cout << boost::math::tgamma(cpp_dec_float_100(3.0/2.0)) << std::endl;

    cpp_dec_float_50 my_pi = boost::math::constants::pi<cpp_dec_float_50>();
    std::cout << "my_pi = " << my_pi << std::endl;
    std::cout << "4atg1 = " << 4*boost::multiprecision::atan(cpp_dec_float_50(1)) << std::endl;
    std::cout << "I1(0) = " << my_pi*my_pi/12 << std::endl;
}

// Вывод значений в файл
void printResultToFile(matrix_type::_vector x, cpp_dec_float_50 k, std::string varName)
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
    //std::setprecision(std::numeric_limits<cpp_dec_float_50>::max_digits10);
    std::ofstream fout; // создаём объект класса ofstream для записи 
    fout.open(fileName);
    fout.precision(std::numeric_limits<cpp_dec_float_50>::max_digits10);
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
                matrix_type::_vector X, cpp_dec_float_50 k)
{
    cpp_dec_float_50 y;
    for (int i = 0; i < I_base.size(); i++) {
        y = log(1 + exp(x0.at(i)));
        I_base.at(i) = pow((I_base.at(i)*exp(x0.at(i)) / y0.at(i)), 1 / k);
    }

    for (int i = 0; i < I_additional.size(); i++) {
        y = log(1 + exp(X.at(i)));
        I_additional.at(i) = pow((I_additional.at(i)*exp(X.at(i)) / Y.at(i)), 1 / k);
    }
}

static void computeIntegral(std::vector<cpp_dec_float_50> x0, 
                            std::vector<cpp_dec_float_50> X,
                            std::vector<cpp_dec_float_50>& I_base, 
                            std::vector<cpp_dec_float_50> &I_additional,
                            cpp_dec_float_50 x_div, int N_gorner, int k)
{
    for (int i = 0; i < x0.size(); i++) {
        if (x0.at(i) > x_div) {
            //cpp_dec_float_50 t = fdsf::get_T_max(X.at(i), k);
            //std::cout << "t = " << t << std::endl;
            //cpp_dec_float_50 t = 60.0;
            //cpp_dec_float_50 t = 75.0;
            cpp_dec_float_50 t = 120.0;
            I_base.push_back(fdsf::Richardson_mesh_refinement(x0.at(i), t, k));
        }
        else {
            I_base.push_back(fdsf::Gorner(x0.at(i), N_gorner, k));
            //I_base.push_back(GornerSchemeForPrecesionY( x0.at(i), N));
        }
        std::cout << "x0: " << x0.at(i) << " I_base: " << I_base.at(i) << std::endl;
    }

    std::cout << "      I1(0) : " << fdsf::PI*fdsf::PI / 12 << std::endl;
    std::cout << "x_div = " << x_div<< std::endl;
    for (int i = 0; i < X.size(); i++) {
        if (X.at(i) > x_div) {
            //cpp_dec_float_50 t = fdsf::get_T_max(X.at(i), k);
            //cpp_dec_float_50 t = 60.0;
            //cpp_dec_float_50 t = 75.0;
            cpp_dec_float_50 t = 120;
            I_additional.push_back(fdsf::Richardson_mesh_refinement(X.at(i), t, k));
        }
        else {
            I_additional.push_back(fdsf::Gorner(X.at(i), N_gorner, k));
        }
        std::cout << "X: " << X.at(i) << " I_add: " << I_additional.at(i) << std::endl;
        //<< std::fixed << std::setprecision(std::numeric_limits<cpp_dec_float_50>::max_digits10)
    }

}
//Получение необходимого числа членов для схемы Горнера произвольной заданной точности
cpp_dec_float_50 get_N_for_Gorner(cpp_dec_float_50 X, cpp_dec_float_50 k)
{
    int i = 1;
    cpp_dec_float_50 Sold = 0, S;
  
    while (true)
    {
        // Итерационное вычисление S
        S = -log(fdsf::epsilon*pow(Sold + 1, k + 1))/abs(X);
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
    const cpp_dec_float_50 one = 1.0;
    matrix_type::_matrix H(N, matrix_type::_vector(N, 0));

    for (size_t n = 0; n < N; n++) {
        for (size_t m = 0; m < N; m++) {
            H[n][m] = one / ((m+1)+(n+1)-1);
        }
    }

    matrix_type::_matrix H_inv = inverse(H);
    std::ofstream fout;
    fout.open("out\\HILBERTTEST.txt");
    for (size_t n = 0; n < N; n++) {
        for (size_t m = 0; m < N; m++) {
            fout << H_inv[n][m] << " ";
        }
        fout << "\n\r";
    }
    fout.close();
}

void testCutCOeffs()
{
    matrix_type::_vector a = { 0.0560148791230902149024568,
        0.0351117957891800867706741,
        0.0021834386943672331415760,
        0.0002464861525522946634693,
        0.0000092228177886669241259 },
        b = { -0.0611726208769112866900252,
        0.0279968542816146833953639,
        -0.0007512148294307540141223,
        0.0000860680747142919882956 };

    std::vector<cpp_dec_float_50> x0, X, y0, Y;
    std::vector<cpp_dec_float_50> I_base, I_additional;
    cpp_dec_float_50 x_div = cpp_dec_float_50(-0.1);
    const int N_gorner = 461; //S = get_N_for_Gorner(x_div, k);
    const int k = 4;
    const int N_base = 4;
    
    // Расчет значения интеграла в базовых узлах
    fdsf::SetLinearTrigonometricGrid(y0, x0, Y, X, N_base);

    matrix_type::_vector delta_base(y0.size(), 0), delta_add(Y.size(), 0);
    // Расчет схемы Горнера и подсчета интеграла на Гауссовой сетке
    computeIntegral(x0, X, I_base, I_additional, x_div, N_gorner, k);
    GetValue_w(I_base, I_additional, y0, x0, Y, X, k);

    GetApproxomateValues(a, b, y0, Y, I_additional, I_base, delta_base, delta_add, N_base);
    printResultToFile(delta_base, k, "delta_base"); printResultToFile(delta_add, k, "delta_add");
}

int main()
{
    std::cout.precision(std::numeric_limits<cpp_dec_float_50>::max_digits10);
    std::vector<cpp_dec_float_50> x0, X, y0, Y;
    std::vector<cpp_dec_float_50> I_base, I_additional;
    cpp_dec_float_50 x_div = cpp_dec_float_50(-0.1);
    std::cout << "x_div = " << x_div << std::endl;
    //int N_gorner = 260, k = 1;
    //int N_gorner = 214, k = 2;
    //int N_gorner = 165, k = 1;
    const int N_gorner = 461; //S = get_N_for_Gorner(x_div, k);
    const int k = 4;
    const int N_base = 4;

    //cpp_dec_float_50 S = get_N_for_Gorner(x_div, k);
    //std::cout << "N_gorner = " << S << std::endl;
    //testOnHilbertMatrix();
//#if 0

    // Расчет значения интеграла в базовых узлах
    fdsf::SetLinearTrigonometricGrid(y0, x0, Y, X, N_base);

    // Расчет схемы Горнера и подсчета интеграла на Гауссовой сетке
    computeIntegral(x0, X, I_base, I_additional, x_div, N_gorner, k);
    //std::cout << "      I1(0)_teor : " << fdsf::PI*fdsf::PI / 12 << std::endl;
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
    // Коэфициенты из статьидля N = 4
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
    //object.fill_matrix(N_base, I_base, y0, B, A);

    //A_inv = inverse(A);
    //object.find_coefficients(A_inv, B, a, b, N_base);
    //printResultToFile(a, k, "a"); printResultToFile(b, k, "b");

    //GetApproxomateValues(a, b, y0, Y, I_additional, I_base, delta_base, delta_add, N_base);
    //printResultToFile(delta_base, k, "delta_base"); printResultToFile(delta_add, k, "delta_add");
//#endif
    testCutCOeffs();
    getchar();
    return 0;
}