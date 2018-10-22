import math
import numpy as np


def trapz(f, a, b, N):
    h = (b - a) / N
    I = (f(a) + f(b)) / 2

    for i in range(1, N):
        I += f(a + i*h)

    return h*I

def func_cos(x):
    return 1.0 / (2 - cos(x))


def func_demo(x):
    n = 0
    a = 0.9
    return cos(n*x)/(1 - 2*a*cos(x) + a*a)


def func_exp_sin(x):
    return exp(sin(x))


def func_exp_cos(x):
    return exp(cos(x))


def Richardson(f, a, b):
    N = 1;
    I = trapz(f, a, b, N)
    a_ = 0.9
    n = 0
    I_prec = (a_**n)*np.pi / (1 - a_*a_ )
    # I_prec = np.pi / ((a_*a_ - 1)*a_)
    print "N = " + N + ": I = " + I 
    std::ofstream fout;
    fout.open("demo.txt");
    do {
        I_2n = trapz(f, a, b, N + 1);

        # stop_criteria = (I / I_2n) - 1;
        stop_criteria = (I_2n / I_prec) - 1;
        I = I_2n;
        N = N + 1;
        # std::cout << "N = " << N << ": I = " << I << std::endl;
        # std::cout << "N = " << N << ": d = " << abs(stop_criteria) << std::endl;
        fout << abs(stop_criteria) << std::endl;
    } while (abs(stop_criteria) > epsilon);

    fout.close()
    return I

def checkTrapz(a, b):
    Richardson(func_demo, a, b)
    # Richardson(func_cos, a, b)
    # Richardson(func_exp_cos, a, b)
