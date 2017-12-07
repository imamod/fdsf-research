#pragma once

#include "BasicTypes.h"
#include "Fdsf.h"
#include <iostream>

/**
 * ������������ ��������� �������
 */
BmpMatrix eye(const size_t N);

// �������� �������� ������� ������� ������
BmpMatrix inverse(const BmpMatrix& A);

/**
* ������ ������� ��� ������������� ������� ����� ��������.
*/
BmpVector solveSystemForIntegerIndex(const BmpVector& z, const BmpVector& y0, size_t N);

/**
* ������� ��������� �����������������x �������� I. �������������� ��� �����
*/
BmpVector approximateFunctionValueIntegerIndex(const BmpVector& a, const BmpVector& b, const BmpVector& y, size_t N_base);

/**
* �������� ����������
*/
BmpVector relativeApproximationError(const BmpVector& precisionF, const BmpVector& approximateF);

/********************************************************
 *                 ��������� �������                    *
 ********************************************************/

/**
 * ������ ������� ��� ������ ������������� ��������� ��������
 */
BmpVector solveRightApproximationSystem(BmpReal k, size_t N, const BmpVector& y0, const BmpVector& I_base);

/**
 * �������� ������ ������������ �������� � ������ � ��� ��������� �������� k, ������ �������������
 */
BmpVector approximateValueRight(const BmpVector& a, const BmpVector& b, const BmpVector& y, BmpReal k);
