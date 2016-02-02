#pragma once

#include <math.h>
#include <vector>

const double PI = 3.141592653589793238463;

// �������� �����
void setLinMesh(std::vector<double> &y0, std::vector<double> &x0, std::vector<double> &Y, std::vector<double> &X, int N_base);

// ����������� �����
void setTrigMesh(std::vector<double> &y0, std::vector<double> &x0, std::vector<double> &Y, std::vector<double> &X, int N_base);

// ��������� �����
void setLinTrigMesh(std::vector<double> &y0, std::vector<double> &x0, std::vector<double> &Y, std::vector<double> &X, int N_base);