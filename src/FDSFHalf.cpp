#include "BasicService.h"
#include "Fdsf.h"

namespace fdsf {

    using FermiFunction = std::function<BmpReal(const BmpReal& ksi, const BmpReal& x, BmpReal k, BmpReal a, const integration_segment_values& isv)>;

    BmpReal fermi_dirak_half_integer(BmpReal ksi, BmpReal x, 
                                      BmpReal k, BmpReal a, 
                                      integration_segment_values isv) {
        BmpReal denom = (1 + ksi) * (BmpReal(isv.N - isv.n) / BmpReal(isv.N));
        //BmpReal denom = (1 - ksi * ksi);
        BmpReal exp_ksi = exp(-a * ksi * ksi / denom);

        return (2*pow(a, k + 1) * pow(ksi, 2 * k + 1) * exp_ksi) /
               (pow(denom, k + 2) * (exp_ksi + exp(-x)));
    }

    BmpReal fermi_dirak_m3half(BmpReal ksi, BmpReal x,
                                BmpReal k, BmpReal a,
                                integration_segment_values isv) {
        BmpReal exp_ksi = exp(-a * ksi * ksi / (1 - ksi * ksi));
        BmpReal sum_exp = exp_ksi + exp(-x);
        BmpReal exp_diff = exp( -a * ksi * ksi / (1 - ksi * ksi) - x);

        return (-4 * exp_diff * sqrt(abs(a)) * exp_ksi) /
            (pow(1 - ksi * ksi, BmpReal( 3.0 / 2 )) * sum_exp * sum_exp);
    }

    // Euler-Macloren Formulas
    static BmpReal trapz(FermiFunction f, BmpReal x, const BmpReal k, size_t N, BmpReal a) {
        BmpReal h = BmpReal(1.0 / N);
        integration_segment_values isv = {0, N};
        BmpReal u0 = f(0, x, k, a, isv);
        // uN принудительно задаем нулем, чтобы не было переполнени€
        BmpReal I = u0 / 2; 
//#if 0
        // true work
        for (size_t i = 1; i < N; i++) {
            isv.n = i;
            I += f(i * h, x, k, a, isv);
        }
//#endif
#if 0
        // TODO: расчет на 2 узлах одновременно??? проверить
        for (size_t i = 1; i < N / 2; i = i + 2) {
            I += f(i*h, x, k, a) + f((N - i)*h, x, k, a);
        }

        if (N == 2) {
            I += f(N*h/2, x, k, a);
        }
#endif
        return h*I;
    }

    static BmpReal quad(FermiFunction f, BmpReal x, const BmpReal k,
                         size_t N, BmpReal a) {
        BmpReal I = 0;
        BmpReal h = BmpReal(1.0 / N);
        for (size_t i = 0; i < N; i++) {
            //I += f((i + 1.0/2)*h, x, k, a);
        }
        I *= h;
        return I;
    }

    static BmpReal simpson(FermiFunction f, const BmpReal& x, const BmpReal& k, size_t N, BmpReal a) {
        BmpReal I = 0;
        BmpReal h = BmpReal(1.0 / N);
        for (size_t i = 0; i < N; i++) {
            //I += f((i + 1.0 / 2)*h, x, k, a);
        }
        I *= h;

        return I;
    }

    BmpReal euler_maclaurin_method(BmpReal x, const BmpReal k, int N, BmpReal& a) {
        //BmpReal a = newton::NewtonsMethod(x, k);
        a = NewtonsMethod(x, k);
        if (k == -3.0 / 2) {
            return trapz(fermi_dirak_m3half, x, k, N, a);
        } 
        else {
            return trapz(fermi_dirak_half_integer, x, k, N, a);
        }
    }

    // TODO: split this function
    void SetLinearTrigonometricGridRight(BmpVector &y_base,
                                         BmpVector &x_base,
                                         BmpVector &Y,
                                         BmpVector &X, size_t N_base) {
        using namespace fdsf;
        size_t n_additional = 11;
        const BmpReal alpha = 2 / (2 + PI);
        const BmpReal one = BmpReal(1);
        const BmpReal num2 = BmpReal(2); //if integer
        const BmpReal x_star = BmpReal(3);
        const BmpReal y_star = BmpReal(log(1 + exp(x_star))); // if half-integer
                                                              //BmpReal baseSize = BmpReal(2 * N_base + 1); // if integer || half-integer & !fixed a(N+1)
        BmpReal baseSize = BmpReal(2 * N_base); // if half-integer & fixed a(N+1)
                                                //BmpReal baseSize = BmpReal(N_base); // if poly approximation

        //const BmpReal y_star_inv = 1 / (y_star * y_star);
        const BmpReal y_star_inv = 1 / y_star;
        //const BmpReal y_star_inv = 1 / pow(y_star, 0.5);
        //const BmpReal y_star_inv = 1 / pow(y_star, 0.25 );
        //const BmpReal y_star_inv = 1 / pow(y_star, 3.0 / 2);

        // «адаютс€ базовые узлы интерпол€ции
        for (size_t j = 1; j <= baseSize; j++) {
            y_base.push_back(y_star_inv / num2*(num2 * alpha*j / baseSize
                + (one - alpha)*(one - cos(PI*j / baseSize))));
        }

        // «адаютс€ дополнительные точки
        Y.push_back(y_base[0] / n_additional);

        for (size_t i = 1; i < n_additional; i++) {
            Y.push_back(Y[i - 1] + y_base[0] / n_additional);
        }

        for (size_t index = 1; index < y_base.size(); index++) {
            for (size_t i = 0; i < n_additional; i++) {
                Y.push_back(Y.back() + (y_base[index] - y_base[index - 1]) / n_additional);
            }
        }

        // –азворачиваем y
        std::reverse(y_base.begin(), y_base.end());
        for (size_t j = 0; j < baseSize; j++) {
            //y_base[j] = 1.0 / pow(y_base[j], 0.5);
            y_base[j] = 1.0 / y_base[j];
            //y_base[j] = 1.0 / (y_base[j] * y_base[j]);
            //y_base[j] = 1.0 / (pow(y_base[j], 4));
            //y_base[j] = 1.0 / (pow(y_base[j], 2.0 / 3));
            x_base.push_back(log(exp(y_base[j]) - one));
        }

        std::reverse(Y.begin(), Y.end());
        //std::cout << Y.size() << std::endl;
        for (size_t j = 0; j < Y.size(); j++) {
            //Y[j] = 1.0 / pow(Y[j], 0.5);
            Y[j] = 1.0 / Y[j]; 
            //Y[j] = 1.0 / (Y[j] * Y[j]);
            //Y[j] = 1.0 / pow(Y[j], 4);
            //Y[j] = 1.0 / pow(Y[j], 2.0 / 3);
            X.push_back(log(exp(Y[j]) - one));
        }
    }

}; //fdsf