#pragma once

#include "BasicTypes.h"

/**
 * Сформировать единичную матрицу
 */
BmpMatrix eye(const size_t N);

// Получить обратную матрицу методом Гаусса
BmpMatrix inverse(const BmpMatrix& A);

/**
* Решить систему для аппроксимации функций целых индексов.
*/
BmpVector solveSystemForIntegerIndex(const BmpVector& z, const BmpVector& y0, size_t N);

/**
* Функция получения аппроксимированныx значений I. Использовалось для целых
*/
BmpVector approximateFunctionValueIntegerIndex(const BmpVector& a, const BmpVector& b, const BmpVector& y, size_t N_base);

/**
* Получить погрешноть
*/
BmpVector relativeApproximationError(const BmpVector& precisionF, const BmpVector& approximateF);

/********************************************************
 *                 Полуцелые индексы                    *
 ********************************************************/

/**
 * Решить систему для правой аппроксимации полуцелых индексов
 */
BmpVector solveRightApproximationSystem(BmpReal k, size_t N, const BmpVector& y0, const BmpVector& I_base);

/**
 * Получить вектор приближенных значений в точках у для полуцелых индексов k, правая аппроксимация
 */
BmpVector approximateValueRight(const BmpVector& a, const BmpVector& b, const BmpVector& y, BmpReal k);
