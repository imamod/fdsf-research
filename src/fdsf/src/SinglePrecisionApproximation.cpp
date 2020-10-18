/**
 * Реализация вычисления функций ФД с помощью формул 2-кусковой аппроксимации
 * с точностью single precision
 */
#include "SinglePrecisionApproximation.h"
#include "FdIndex.h"
#include "Constants.h"
#include "Gamma.h"
#include "Errors.h"

namespace {

    // Корректное вычисление функции ФД нулевого индекса во всем диапазоне аргумета x
    BmpReal fd_0(BmpReal x) {
        BmpReal exp_x = exp(x);
        if (exp_x < 1e-8) {
            return BmpReal(2 * exp_x) / (2 + exp_x);
        } else if (exp_x < 1e+8) {
            return log(1 + exp_x);
        }
        return x + BmpReal(2) / (1 + 2 * exp_x);
    }

    // Возвращает коэффициенты числителя левой аппроксимации для заданного индекса
    BmpVector leftNomParamsByIndex(BmpReal k) {
        if (fdsf::index::M3_HALF == k) {
            return { 1,

            };
        } else if (fdsf::index::M1_HALF == k) {
            return { 1,
                     0.49870870194718009,
                     0.11316839723076555,
                     0.014666169008251018,
                     0.001032378227563413
            };
        } else if (fdsf::index::P1_HALF == k) {
            return { 1,
                     0.3820834770485817,
                     0.080007308526546694,
                     0.010045334815231399,
                     0.00067271857491846276
            };
        } else if (fdsf::index::P3_HALF == k) {
            return { 1,
                     0.3061333832629316,
                     0.065014952751880628,
                     0.0079062467198127706,
                     0.00050989014982860681
            };
        } else if (fdsf::index::P5_HALF == k) {
            return { 1,
                     0.25572368082794128,
                     0.056222494510620891,
                     0.0067854625228846999,
                     0.00044476579091679014
            };
        } else if (fdsf::index::P7_HALF == k) {
            return { 1,
                     0.20556389789635432,
                     0.046771068015914352,
                     0.0054384771639206519,
                     0.00036630399605996899
            };
        }
        throw UnsupportedFdFunction();
    }

    // Возвращает коэффициенты знаменателя левой аппроксимации для заданного индекса
    BmpVector leftDenomParamsByIndex(BmpReal k) {
        if (fdsf::index::M3_HALF == k) {
            return { 1,
                // TODO:
            };
        } else if (fdsf::index::M1_HALF == k) {
            return { 1,
                     0.084497999356244691,
                     0.023296798173987554,
                     0.0011595586921657741
            };
        } else if (fdsf::index::P1_HALF == k) {
            return { 1,
                     0.089191759958339389,
                     0.021303987930878066,
                     0.0010768539618197792
            };
        } else if (fdsf::index::P3_HALF == k) {
            return { 1,
                     0.090652313530881656,
                     0.021057720402495761,
                     0.0010291094073409113
            };
        } else if (fdsf::index::P5_HALF == k) {
            return { 1,
                     0.091080083402630407,
                     0.021688411740569791,
                     0.0010773768482295054
            };
        } else if (fdsf::index::P7_HALF == k) {
            return { 1,
                     0.075334948694944615,
                     0.021126068076682714,
                     0.0010018164617804359
            };
        }
        throw UnsupportedFdFunction();
    }

    // Возвращает коэффициенты числителя правой аппроксимации для заданного индекса
    BmpVector rightNomParamsByIndex(BmpReal k) {
        if (fdsf::index::M3_HALF == k) {
            return { 1,
                     48.598302155831334,
                     5442.2885034815845,
                     291822.15014591813,
                    -764104.49450398982
            };
        } else if (fdsf::index::M1_HALF == k) {
            return { 1,
                     29.502444655512591,
                     4185.4823745111953,
                     96683.277205527993,
                     1831673.6164128184
            };
        } else if (fdsf::index::P1_HALF == k) {
            return { 1,
                     64.829801157931797,
                     4616.5501019326621,
                     98092.21572034806,
                     389335.91363818944
            };
        } else if (fdsf::index::P3_HALF == k) {
            return { 1,
                     79.962518648247169,
                     5068.0728371160803,
                     115844.87551961094,
                     372013.5346448943
            };
        } else if (fdsf::index::P5_HALF == k) {
            return { 1,
                     86.738516205206906,
                     5207.1771901756874,
                     132850.77403663844,
                     439660.24346425384
            };
        } else if (fdsf::index::P7_HALF == k) {
            return { 1,
                     109.68251185932604,
                     5993.7024718266912,
                     201633.54644721746,
                     729559.12439405918
            };
        }
        throw UnsupportedFdFunction();
    }

    // Возвращает коэффициенты знаменателя правой аппроксимации для заданного индекса
    BmpVector rightDenomParamsByIndex(BmpReal k) {
        if (fdsf::index::M3_HALF == k) {
            return { 1,
                     39.062015561877615,
                     5373.2838163729757,
                     109639.78241016809,
                     1755046.2371845618
            };
        } else if (fdsf::index::M1_HALF == k) {
            return { 1,
                     24.048332989117625,
                     4041.1997121249742,
                     59618.305730951251,
                     1066504.1470012292
            };
        } else if (fdsf::index::P1_HALF == k) {
            return { 1,
                     62.094313335226616,
                     4451.6669005131698,
                     80941.512053528801,
                     422893.86418903247
            };
        } else if (fdsf::index::P3_HALF == k) {
            return { 1,
                     79.214104078915625,
                     5019.081662083976,
                     110655.69410856254,
                     472970.26044385135
            };
        } else if (fdsf::index::P5_HALF == k) {
            return { 1,
                     87.309336705484384,
                     5254.8372477423982,
                     134339.21008465067,
                     578488.72409081459
            };
        } else if (fdsf::index::P7_HALF == k) {
            return { 1,
                     110.91242147890443,
                     6114.0387978623621,
                     206573.59537407756,
                     952009.96586287022
            };
        }
        throw UnsupportedFdFunction();
    }
}

// Левая аппроксимация
BmpReal single_precision_formula::left(BmpReal k, BmpReal x) {
    BmpReal y = fd_0(x);
    BmpVector params = leftNomParamsByIndex(k);
    BmpReal nominator = 0;
    for (size_t i = 0; i < params.size(); ++i) {
        nominator += params.at(i) * pow(y, i);
    }
    params = leftDenomParamsByIndex(k);
    BmpReal denominator = 0;
    for (size_t i = 0; i < params.size(); ++i) {
        denominator += params.at(i) * pow(y, i);
    }
    return factorial(k) * y * pow(nominator / denominator, k);
}

// Правая аппроксимация
BmpReal single_precision_formula::right(BmpReal k, BmpReal x) {
    BmpReal y = fd_0(x);
    BmpVector params = rightNomParamsByIndex(k);
    BmpReal nominator = 0;
    for (size_t i = 0; i < params.size(); ++i) {
        nominator += params.at(i) * pow(y, -2*i);
    }
    params = rightDenomParamsByIndex(k);
    BmpReal denominator = 0;
    for (size_t i = 0; i < params.size(); ++i) {
        denominator += params.at(i) * pow(y, -2*i);
    }
    const BmpReal C2 = (k + 1) * pow(pi(), 2) / 3;
    return pow(nominator / denominator * C2 + pow(y, 2), k / 2) * y / (k + 1);
}
