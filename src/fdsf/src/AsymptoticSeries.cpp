#include "AsymptoticSeries.h"
#include "Constants.h"
#include "FdIndex.h"
#include "Logger.h"

#include <iostream>

namespace asympt_series {

    /* Коэффициенты асимптотического рядадля ФД полуцелого индекса */

    const BmpVector A_m3half = {
        1.2337005501361697, 12.42980588715135, 320.15011231365116, 15776.478129534085, 1277607.759031492, 154373786.61820173, 26055305431.853966,
        5856196906551.461, 1.6909961741967747e+15, 6.10028613757551e+17, 2.6887030375639702e+20, 1.421651985282955e+23
    };

    const BmpVector A_mhalf = {
        -0.41123351671205655, -1.7756865553073355, -29.104555664877378, -1051.7652086356056, -67242.51363323643, -6711903.766008771, -965011312.2908876,
        -188909577630.6923, -48314176405622.13, -1.564175932711669e+16, -6.252797761776675e+18, -3.0247914580488403e+21
    };

    const BmpVector A_half = {
        1.2337005501361697, 1.0654119331844014, 9.701518554959126, 242.7150481466782, 11866.325935277016, 958843.3951441102, 115801357.4749065,
        19542370099.726788, 4392197855056.5576, 1.2682507562527045e+15, 4.5752178744707386e+17, 2.0165276386992272e+20
    };

    const BmpVector A_3half = {
        6.168502750680848, -1.7756865553073355, -6.92965611068509, -110.32502188485374, -3955.4419784256725, -252327.20924845006, -25174208.146718808,
        -3618957425.875331, -708419008880.09, -181178679464672.06, -5.865663941629151e+16, -2.3447995798828216e+19
    };

    const BmpVector A_5half = {
        14.39317308492198, 12.42980588715135, 9.701518554959126, 85.80835035488624, 2129.8533729984388, 103899.43910230296, 8391402.715572935,
        1013308079.2450927, 170997691798.64243, 38431841098566.81, 1.1097202051730828e+16, 4.0033163558975007e+18
    };

    const BmpVector A_7half = {
        25.907711552859563, 111.86825298436214, -29.104555664877378, -110.32502188485374, -1742.607305180541, -62339.66346138178, -3974874.970534548,
        -396511857.09590584, -56999230599.54747, -11157631286680.688, -2.8535662418736415e+15, -9.238422359763464e+17
    };

    // Получение x_min и N_max для конкретного индекса
    Limits limits(BmpReal k) {
        if (fdsf::index::M3_HALF == k) {
            return{ 52, 11 };
        } else if (fdsf::index::M1_HALF == k) {
            return{ 39, 10 };
        } else if (fdsf::index::P1_HALF == k) {
            return{ 35, 10 };
        } else if (fdsf::index::P3_HALF == k) {
            return{ 33, 10 };
        } else if (fdsf::index::P5_HALF == k) {
            return{ 30, 10 };
        } else if (fdsf::index::P7_HALF == k) {
            return{ 29, 10 };
        }
        throw std::invalid_argument("Unsuppported k-index value");
    }

    // Получает коэффициенты асимптотического ряда для конкретного k
    BmpVector coefficents(BmpReal k) {
        if (fdsf::index::M3_HALF == k) {
            return A_m3half;
        } else if (fdsf::index::M1_HALF == k) {
            return A_mhalf;
        } else if (fdsf::index::P1_HALF == k) {
            return A_half;
        } else if (fdsf::index::P3_HALF == k) {
            return A_3half;
        } else if (fdsf::index::P5_HALF == k) {
            return A_5half;
        } else if (fdsf::index::P7_HALF == k) {
            return A_7half;
        }
        throw std::invalid_argument("Unsupported k-index value");
    }

    // Вычислить значение ФД для x >= x_min
    // Схема Горнера для асимптотического ряда (см. формулу (32) препринт 2)
    BmpReal calculate(BmpReal k, BmpReal x) {
        Logger log("asymp::calculate");
        // Получаем коэффициенты для конкретного k
        BmpVector A = coefficents(k);
        // Получаем предельные значения для каждого k
        Limits data = limits(k);
        BmpReal x_2_m1 = 1.0 / (x*x);
        BmpReal sum = A.back() * x_2_m1;
        // TODO: когда поймем, какое N_max, использовать его
        //A.size() - 2
        for (int n = data.N_max; n > -1; --n) {
            sum = (A.at(n) + sum) * x_2_m1;
            if (n > 3) {
                //std::cout << "A(" << n + 2 << ")/A(" << n + 1 << ") = "
                   // << (A.at(n + 1) / A.at(n))*x_2_m1 << std::endl;
            }
        }
        // Главный член асимптотики
        BmpReal mainPart = pow(x, k + 1) / (k + 1);
        return mainPart * (sum + 1);
    }
}
