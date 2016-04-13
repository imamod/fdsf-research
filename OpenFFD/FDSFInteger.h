#pragma once

#include <vector>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/detail/default_ops.hpp>
#include <boost/multiprecision/number.hpp>
#include <boost/math/special_functions/gamma.hpp>

namespace fdsf{
    using namespace boost::multiprecision;

    // Желаемая точность расчета
    const cpp_dec_float_50 epsilon = 1e-25;
    
    // Значение pi
    const cpp_dec_float_50 PI = boost::math::constants::pi<cpp_dec_float_50>();

    // Значение Ik(0) для индекса k = 0,1,2,3
    const static cpp_dec_float_50 I_k_0[] = { log(cpp_dec_float_50(2)),
                                              PI*PI / 12.0, 
                                              1.8030853547393952,// значение из статьи
                                              7.0*PI*PI*PI*PI / 120.0 };

    // Задает линейно-тригонометрическое сетку в базовых узлах, плюс 10(?)
    // дополнительных точек между каждой парой базовых узлов.
    void SetLinearTrigonometricGrid(std::vector<cpp_dec_float_50> &y_base, 
                                    std::vector<cpp_dec_float_50> &x_base, 
                                    std::vector<cpp_dec_float_50> &Y, 
                                    std::vector<cpp_dec_float_50> &X, int N_base);

    // Вычисление Г-функции
    // TODO: сделать для полуцелых индексов
    // TODO: лучше подрубить boost
    cpp_dec_float_50 factorial(cpp_dec_float_50 k);

    // Вычисляет значение функции ФД индекса k = 1, 2, 3 в точке x при заданном
    // значении t. Функция представлена в виде t^k /(exp(x)+exp(t))
    cpp_dec_float_50 FermiDirakFunction(cpp_dec_float_50 t, cpp_dec_float_50 x, cpp_dec_float_50 k);

    // Вычисляет значение функции ФД индекса k = 1, 2, 3 в точке x по схеме 
    // Горнера при x <= -0.1. N - число членов в схеме Горнера для достижения 
    // машинной точности
    cpp_dec_float_50 Gorner(cpp_dec_float_50 x, int N, cpp_dec_float_50 k);

    // Алгоритм выбора T_max длякаждого x (а не фиксированные для всех)
    // NB: для новой формы записи подынтегральной функции
    cpp_dec_float_50 get_T_max(cpp_dec_float_50 X, int k);

    // Вычисляет значение функции ФД индекса k = 1, 2, 3 в точке x по формуле 
    // Гаусса-Кристоффеля с 5 узлами на равномерной сетке. 
    // По статье, для достижения требуемой точности на отрезке -0.1<x<0 нужно
    // брать для k = 0 : Tmax = 50, 
    //           k = 1 : Tmax = 60,
    //           k = 2 : Tmax = 75, 
    //           k = 3 : Tmax = 100. 
    // TODO: std::function
    cpp_dec_float_50 FDGK5(cpp_dec_float_50(*Ft)(cpp_dec_float_50, cpp_dec_float_50, cpp_dec_float_50), 
                 cpp_dec_float_50 x, cpp_dec_float_50 T, cpp_dec_float_50 k, int N);

    // Сгущение по Ричардсону результата работы функции FDGK5
    cpp_dec_float_50 Richardson_mesh_refinement(cpp_dec_float_50 x, cpp_dec_float_50 t, cpp_dec_float_50 k);
    
    namespace integer {
        // Вычисляет значение функции ФД индекса k=1 в точке x
        cpp_dec_float_50 FD_I1(cpp_dec_float_50 x);

        // Вычисляет значение функции ФД индекса k=2 в точке x
        cpp_dec_float_50 FD_I2(cpp_dec_float_50 x);

        // Вычисляет значение функции ФД индекса k=3 в точке x
        cpp_dec_float_50 FD_I3(cpp_dec_float_50 x);
    } // namespace integer

} // namespace fdsf