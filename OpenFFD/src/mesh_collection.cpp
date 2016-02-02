// Implementation of difference_scheme.h
//
#include "stdafx.h"
#include "mesh_collection.h"

void setLinMesh(std::vector<double> &y0, std::vector<double> &x0, std::vector<double> &Y, std::vector<double> &X, int N_base)
{
    int n_additional = 11;

    // Задается линейная сетка
    for (int j = 1; j <= 2 * N_base + 1; j++) {
        y0.push_back(j*log(2) / (2 * N_base + 1));
        x0.push_back(log(exp(y0.at(j - 1)) - 1));
    }

    Y.push_back(y0.at(0) / n_additional);
    X.push_back(log(exp(Y.at(0)) - 1));

    for (int i = 1; i < n_additional; i++)
    {
        Y.push_back(Y.at(i - 1) + y0.at(0) / n_additional);
        X.push_back(log(exp(Y.at(i)) - 1));
    }

    for (int index = 1; index < y0.size(); index++) {
        for (int i = 0; i < n_additional; i++) {
            Y.push_back(Y.back() + (y0.at(index) - y0.at(index - 1)) / n_additional);
            X.push_back(log(exp(Y.back()) - 1));
        }
    }

}

void setTrigMesh(std::vector<double> &y0, std::vector<double> &x0, std::vector<double> &Y, std::vector<double> &X, int N_base)
{
    int n_additional = 11;

    // Задается тригонометрическая сетка
    for (int j = 1; j <= 2 * N_base + 1; j++) {
        y0.push_back(0.5*log(2)*(1 - cos(PI*j / (2 * N_base + 1))));
        x0.push_back(log(exp(y0.at(j - 1)) - 1));
    }

    Y.push_back(y0.at(0) / n_additional);
    X.push_back(log(exp(Y.at(0)) - 1));

    for (int i = 1; i < n_additional; i++)
    {
        Y.push_back(Y.at(i - 1) + y0.at(0) / n_additional);
        X.push_back(log(exp(Y.at(i)) - 1));
    }

    for (int index = 1; index < y0.size(); index++) {
        for (int i = 0; i < n_additional; i++) {
            Y.push_back(Y.back() + (y0.at(index) - y0.at(index - 1)) / n_additional);
            X.push_back(log(exp(Y.back()) - 1));
        }
    }

}

void setLinTrigMesh(std::vector<double> &y0, std::vector<double> &x0, std::vector<double> &Y, std::vector<double> &X, int N_base)
{
    int n_additional = 11;
    double alpha = 2 / (2 + PI);

    // Задается линейно-тригонометрическая сетка
    for (int j = 1; j <= 2 * N_base + 1; j++) {
        y0.push_back(0.5*log(2)*(2 * alpha*j / (2 * N_base + 1) + (1 - alpha)*(1 - cos(PI*j / (2 * N_base + 1)))));
        x0.push_back(log(exp(y0.at(j - 1)) - 1));
    }

    Y.push_back(y0.at(0) / n_additional);
    X.push_back(log(exp(Y.at(0)) - 1));

    for (int i = 1; i < n_additional; i++)
    {
        Y.push_back(Y.at(i - 1) + y0.at(0) / n_additional);
        X.push_back(log(exp(Y.at(i)) - 1));
    }

    for (int index = 1; index < y0.size(); index++) {
        for (int i = 0; i < n_additional; i++) {
            Y.push_back(Y.back() + (y0.at(index) - y0.at(index - 1)) / n_additional);
            X.push_back(log(exp(Y.back()) - 1));
        }
    }

}
