#include "Common.h"
#include "Grid.h"
#include "MatrixUtils.h"
#include "Fdsf.h"
#include <algorithm>

namespace {

    // Вспомогательная структура с информацией о узле интерполяции
    struct NodeInfo {
        BmpReal value;    // значение точки
        bool isBaseNode;  // является ли базовым
    };

    BmpVector convNodeInfoToVec(const std::vector<NodeInfo>& v) {
        BmpVector result;
        for (auto it : v) {
            result.emplace_back(it.value);
        }
        return result;
    }

    // Начальное распределение точек
    BmpVector initGrid(size_t N) {
        BmpVector grid;
        for (size_t n = 0; n < N; ++n) {
            grid.emplace_back((2 * n - N) / N);
        }
        return grid;
    }

    // Получает локальный экстремум на заданном промежутке
    template <typename It>
    BmpReal getLocalExtremum(It first, It last, bool isPositive) {
        return isPositive ? *std::max_element(first, last) : *std::min_element(first, last);
    }

    const size_t ADDITIONAL_DOTS_COUNT = 20;

    /**
     * Возвращает массив локальных экстремумов между базовыми узлами интерполяции
     * @param - профиль погрешности
     * @param - число базовых узлов интерполяции
     * @return - массив локальных экстремумов
     */
    BmpVector localInterpolationExtremums(const BmpVector& errorProfile, size_t N) {
        // Рачет экстремумов
        size_t middleDot = ADDITIONAL_DOTS_COUNT / 2;
        BmpVector localExtremums;
        BmpReal firstExtremum = getLocalExtremum<BmpVector::const_iterator>(errorProfile.begin(), errorProfile.begin() + ADDITIONAL_DOTS_COUNT, errorProfile[middleDot] > 0);
        localExtremums.emplace_back(firstExtremum);
        for (size_t i = 1; i < N; ++i) {
            size_t span = i * ADDITIONAL_DOTS_COUNT;
            auto startPos = errorProfile.begin() + span;
            auto endPos = errorProfile.begin() + span + ADDITIONAL_DOTS_COUNT;
            BmpReal currentExtremum = getLocalExtremum<BmpVector::const_iterator>(startPos, endPos, errorProfile[span + middleDot] > 0);
            localExtremums.emplace_back(currentExtremum);
        }
        return localExtremums;
    }

    bool isBaseNode(const NodeInfo& x) {
        return x.isBaseNode;
    }

    BmpVector localInterpolationExtremumsWithIndependentAddDots(const std::vector<NodeInfo>& errorProfile, size_t N) {
        BmpVector localExtremums;
        auto startPos = errorProfile.begin() + 1;
        auto endPos = std::find_if(startPos, errorProfile.end(), isBaseNode);
        auto middleDot = std::distance(errorProfile.begin(), endPos) / 2;
        BmpVector localPart(convNodeInfoToVec(std::vector<NodeInfo>{errorProfile.begin(), endPos}));
        BmpReal firstExtremum = getLocalExtremum<BmpVector::const_iterator>(localPart.begin(), localPart.end(), errorProfile[middleDot].value > 0);
        localExtremums.emplace_back(firstExtremum);
        localPart.clear();
        for (size_t i = 1; i < N; ++i) {
            startPos = endPos;
            endPos = std::find_if(startPos + 1, errorProfile.end(), isBaseNode);
            BmpVector localPart(convNodeInfoToVec(std::vector<NodeInfo>{startPos, endPos}));
            middleDot = std::distance(errorProfile.begin(), startPos) + std::distance(startPos, endPos) / 2;
            BmpReal currentExtremum = getLocalExtremum<BmpVector::const_iterator>(localPart.begin(), localPart.end(), errorProfile[middleDot].value > 0);
            localExtremums.emplace_back(currentExtremum);
        }
        return localExtremums;
    }

    /*class MeltipliersHaveTheSameSign : public std::exception {
        public:
            const char* what() {

            }
    };*/

    /**
     * Возвращает массив скоростей базовых узлов
     * @param  вектор узлов
     * @param  вектор экстремумов
     * @return массив скоростей базовых узлов
     */
    BmpVector getBaseNodesVelocity(const BmpVector& x, const BmpVector& extremums) {
        BmpVector v = { 0 };
        // настроечная константа скорости, рабочий вариант
        const BmpReal MUL_CONST = 1 / (2 * sqrt(3));
        for (size_t i = 0; i < extremums.size() - 1; ++i) {
            // если в списке экстремумов есть соседи одного знака, срыв алгоритма
            if (extremums[i] * extremums[i + 1] > 0) {
                throw std::invalid_argument("Multipliers have the same sign");
            }
            // добавляем скорости движения узлов
            BmpReal v_pm = (extremums.at(i + 1) + extremums.at(i)) / (extremums.at(i + 1) - extremums.at(i));
            v.emplace_back((x.at(i + 2) - x.at(i)) * v_pm * MUL_CONST);
        }
        return v;
    }

    /**
     * Возвращает шаг tau для текущей сетки
     */
    BmpReal getTau(const BmpVector& x, const BmpVector& v) {
        // Условно большое tau
        const BmpReal TAU_BIG = 1e10;
        // расчет tau для каждого узла
        BmpVector tau;
        for (size_t i = 0; i < v.size() - 1; ++i) {
            BmpReal denom = v.at(i) - v.at(i + 1);
            BmpReal nextTau = (denom < 0) ? TAU_BIG : (x.at(i + 1) - x.at(i)) / denom;
            tau.emplace_back(nextTau);
        }
        // настроечный коэффициент tau
        const BmpReal A = 0.2;
        BmpReal tauOptimal = *std::min_element(tau.begin(), tau.end()) * A;
        const size_t MAX_TAU_VALUE = 1;
        if (tauOptimal > MAX_TAU_VALUE) {
            tauOptimal = MAX_TAU_VALUE;
        }
        return tauOptimal;
    }

    /**
     * Сдвиг сетки
     * @param  вектор базовых узлов
     * @param  вектор экстремумов
     * @return вектор сдвинутой сетки
     */
    BmpVector shiftBaseNodesGrid(const BmpVector& x, const BmpVector& extremums) {
        // Скорости движения соседних узлов
        BmpVector v = getBaseNodesVelocity(x, extremums);
        BmpReal tau = getTau(x, v);
        BmpVector shiftedBaseNodes;
        for (size_t i = 0; i < x.size() - 1; ++i) {
            shiftedBaseNodes.emplace_back(x.at(i) + v.at(i) * tau);
        }
        // Последний узел отрезка неподвижен
        shiftedBaseNodes.emplace_back(x.back());
        return shiftedBaseNodes;
    }

    /**
     * param сетка базовых узлов
     * param сетка дополнительных узлов
     * return аппроксимированное значение функции в узлах
     */
    BmpVector approximatePolinom(const BmpVector& x, const BmpVector& X) {

        const size_t baseSize = x.size();
        // Аппроксимируемая функция = 0
        BmpVector U_prec(baseSize, 0);
        // Правая часть
        BmpVector rightPart;
        for (size_t i = 0; i < x.size(); ++i) {
            rightPart.emplace_back(U_prec.at(i) - pow(x.at(i), baseSize));
        }
        // Левая часть
        BmpMatrix A(baseSize, BmpVector(baseSize, 0));
        for (size_t i = 0; i < baseSize; ++i) {
            for (size_t j = 0; j < baseSize; ++j) {
                A.at(i).at(j) = pow(x.at(i), j);
            }
        }
        BmpMatrix A_inv = inverse(A);
        // Решаем СЛАУ A*E = B
        BmpVector E(baseSize, 0);
        for (size_t i = 0; i < baseSize; ++i) {
            for (size_t j = 0; j < baseSize; ++j) {
                E[i] += A_inv[i][j] * rightPart[j];
            }
        }
        // Ограничения на многочлен
       // E.emplace_back(1);
        BmpVector U;
        for (size_t i = 0; i < X.size(); ++i) {
            BmpReal u = pow(X.at(i), baseSize);
            for (size_t j = 0; j < baseSize; ++j) {
                u = u + E.at(j) * pow(X.at(i), j);
            }
            U.emplace_back(u);
        }
        return U;
    }

    // Равномерное распределение точек размера N + 1
    // N - число экстремумов
    BmpVector getUniformDistributionDots(BmpReal N) {
        BmpVector x;
        for (size_t n = 0; n <= N; ++n) {
            x.push_back(BmpReal(2 * n - N) / N);
        }
        return x;
    }

    // Литр справа, Т - число экстремумов
    BmpVector getLinTrigDistribution(BmpReal N) {
        Grid grid(N + 1, ADDITIONAL_DOTS_COUNT);
        grid.setLinearTrigonometricGridRight();
        return grid.base();
    }

    BmpVector powVector(const BmpVector& v, BmpReal p) {
        BmpVector result;
        for (auto it : v) {
            result.emplace_back(pow(it, p));
        }
        return result;
    }

    BmpVector getPrecisionValue(const BmpVector& x) {
        BmpVector result;
        for (auto it : x) {
            result.emplace_back(fdsf::fd_1half(it));
        }
        return result;
    }

    std::vector<NodeInfo> getPrecisionValueWithInfo(const std::vector<NodeInfo>& x) {
        std::vector<NodeInfo> result;
        for (auto it : x) {
            result.emplace_back(NodeInfo{ fdsf::fd_1half(it.value), it.isBaseNode });
        }
        std::cout << "result size " << result.size() << std::endl;
        return result;
    }

    BmpVector getApproximationFormByFDValue(const BmpVector& y, const BmpVector& v) {
        const BmpReal k = 0.5;
        const BmpReal C1 = (k + 1) * k * pow(pi(), 2) / 6;
        BmpVector result;
        for (size_t i = 0; i < v.size(); ++i) {
            result.emplace_back(((v[i] * (k + 1) / pow(y[i], k + 1) - 1) * pow(y[i], 2)) / C1);
        }
        return result;
    }

    BmpVector getFDValueFromApproximation(const BmpVector& Y, const BmpVector& v) {
        BmpVector U;
        const BmpReal k = 0.5;
        const BmpReal C1 = (k + 1) * k * pow(pi(), 2) / 6;
        for (size_t i = 0; i < v.size(); ++i) {
            U.emplace_back((v.at(i) * C1 / pow(Y.at(i), 2) + 1) * (pow(Y.at(i), k + 1) / (k + 1)));
        }
        return U;
    }

    BmpReal calculateI0(BmpReal x) {
        BmpReal exp_x = exp(x);
        if (exp_x < 1e-8) {
            return BmpReal(2 * exp_x) / (2 + exp_x);
        } else if (exp_x < 1e+8) {
            return log(1 + exp_x);
        } else {
            return x + BmpReal(2) / (1 + 2 * exp_x);
        }
    }

    BmpVector getI0Values(const BmpVector& x) {
        BmpVector y;
        for (auto it : x) {
            y.emplace_back(calculateI0(it));
        }
        return y;
    }

    BmpVector approximateRationalFd(const BmpVector& x, const BmpVector& X, const BmpVector& U_additional) {

        BmpVector U_base = getPrecisionValue(x);
        const BmpReal k = 0.5;
        const size_t baseSize = x.size()-1; // т. БЕсконечность в аппроксимации не участвует
        const size_t N = baseSize / 2;
        BmpVector y0 = getI0Values(x);
        BmpVector Y = getI0Values(X);
        const BmpReal C1 = (k + 1) * k * pow(pi(), 2) / 6;
        BmpVector z = getApproximationFormByFDValue(y0, U_base);
        BmpVector I = getApproximationFormByFDValue(Y, U_additional);
        // Правая часть
        BmpVector rightPart;
        for (auto it : z) {
            rightPart.emplace_back(it - 1);
        }
        // Левая часть
        BmpMatrix A(baseSize, BmpVector(baseSize, 0));
        for (size_t i = 0; i < baseSize; ++i) {
            for ( size_t j = 0; j < baseSize; ++j) {
                if (j < N) {
                    A.at(i).at(j) = pow(y0.at(i), -2 * BmpReal(j + 1));
                } else if (j < baseSize) {
                    A.at(i).at(j) = -z.at(i) * pow(y0.at(i), -2 * BmpReal(j + 1 - N));
                }
            }
        }
        BmpMatrix A_inv = inverse(A);
        // Решаем СЛАУ A*E = B
        BmpVector E(baseSize, 0);
        for (size_t i = 0; i < baseSize; ++i) {
            for (size_t j = 0; j < baseSize; ++j) {
                E[i] += A_inv[i][j] * rightPart[j];
            }
        }

        // Раскидываем коэффициенты
        BmpVector a, b;
        for (size_t j = 0; j < E.size(); ++j) {
            if (j < N) {
                a.emplace_back(E.at(j));
            } else if (j <= baseSize) {
                b.emplace_back(E.at(j));
            }
        }

        //
        BmpVector F_additional;
        for (auto y : Y) {
            BmpReal nom = 0, denom = 0;
            for (size_t n = 1; n <= N; ++n) {
                nom += a.at(n-1)*pow(y, (-2 * BmpReal(n)));
            }
            for (size_t m = 1; m <= N; ++m) {
                denom += b.at(m - 1) * pow(y, (-2 * BmpReal(m)));
            }
            F_additional.emplace_back(BmpReal(1.0 + nom) / (1.0 + denom));
        }

        return getFDValueFromApproximation(Y, F_additional);
    }

    // Дополнительные узлы
    // param - вектор базовых узлов
    BmpVector getAdditionalDots(const BmpVector& x) {
        BmpVector additionalDots{ x.front() };
        for (size_t i = 0; i < x.size() - 1; ++i) {
            BmpReal delta = (x.at(i + 1) - x.at(i)) / ADDITIONAL_DOTS_COUNT;
            for (size_t j = 1; j < ADDITIONAL_DOTS_COUNT; ++j) {
                additionalDots.emplace_back(x.at(i) + j * delta);
            }
            additionalDots.emplace_back(x.at(i + 1));
        }
        return additionalDots;
    }

    // Степень величины
    const BmpReal P(2);

    // TODO: рефакторинг
    // Для литр, доп точки нужно получать точки по обратной величине
    std::vector<NodeInfo> getAdditionalDotsWithInfo(const BmpVector& x) {
        BmpVector inverseValue;
        for (const auto& it : x) {
            inverseValue.emplace_back(BmpReal(1) / pow( it, P));
        }
        std::cout << "x size " << x.size() << std::endl;
        std::reverse(inverseValue.begin(), inverseValue.end());
        // Полубесконечный промежуток
        std::vector<NodeInfo> additionalDots{ NodeInfo{x.back(), true} }; // т. бесконечность
        for (size_t i = 0; i < inverseValue.size() - 1; ++i) {
            BmpReal delta = (inverseValue.at(i + 1) - inverseValue.at(i)) / ADDITIONAL_DOTS_COUNT;
            for (size_t j = 1; j < ADDITIONAL_DOTS_COUNT; ++j) {
                additionalDots.emplace_back(NodeInfo{ BmpReal(1) / pow(inverseValue.at(i) + j * delta, 1/ P), false });
            }
            additionalDots.emplace_back(NodeInfo{ BmpReal(1) / pow(inverseValue.at(i + 1), 1 / P), true });
        }
        std::cout << "additionalDots size " << additionalDots.size() << std::endl;
        std::reverse(additionalDots.begin(), additionalDots.end());
        /*
        std::vector<NodeInfo> additionalDots{ { x.front(), true } };
        for (size_t i = 0; i < x.size() - 1; ++i) {
            BmpReal delta = (x.at(i + 1) - x.at(i)) / ADDITIONAL_DOTS_COUNT;
            for (size_t j = 1; j < ADDITIONAL_DOTS_COUNT; ++j) {
                additionalDots.emplace_back(NodeInfo{ x.at(i) + j * delta, false });
            }
            additionalDots.emplace_back(NodeInfo{ x.at(i + 1), true });
        }
        */
        return additionalDots;
    }

    std::vector<NodeInfo> getErrorProfile(const std::vector<NodeInfo>& U_prec, const BmpVector& U) {
        std::vector<NodeInfo> delta;
        for (size_t i = 0; i < U_prec.size(); ++i) {
            delta.emplace_back(NodeInfo{ U[i] / U_prec[i].value - 1, U_prec[i].isBaseNode });
        }
        return delta;
    }

    BmpVector abs(const BmpVector& v) {
        BmpVector absoluteValue;
        for (auto it : v) {
            absoluteValue.push_back(std::abs(it));
        }
        return absoluteValue;
    }

    void testZeroFunction() {
        const BmpReal N = 3;
        BmpVector x = getUniformDistributionDots(N);
        // Значение константы отношения u_max / u_min - 1
        const BmpReal U_RELATION_CONST = 0.01; // Для модельных примеров
        size_t it = 1;
        const size_t MAX_ITERATION_COUNT = 50;
        while (it < MAX_ITERATION_COUNT) {
            std::cout << std::endl << std::endl;
            std::cout << "Current iteration:" << it << std::endl;
            std::cout << "Current grid x: " << std::endl;
            print::vector(x); std::cout << std::endl << std::endl;
            BmpVector additionalDots = getAdditionalDots(x);
            // Аппроксимирующий многочлен (для U = 0 погрешности - сами значения аппроксимации) 
            BmpVector U = approximatePolinom(x, additionalDots);
            BmpVector delta = U;
            // Получаем экстремумы
            BmpVector localExtremums = localInterpolationExtremums(delta, N);
            std::cout << "Current extremums: " << std::endl;
            print::vector(localExtremums); std::cout << std::endl;
            BmpVector moduleLocalExtremums = abs(localExtremums);
            // Диагностика
            BmpReal L = *std::max_element(moduleLocalExtremums.begin(), moduleLocalExtremums.end()) / *std::min_element(moduleLocalExtremums.begin(), moduleLocalExtremums.end()) - 1;
            std::cout << "Current L: " << L << std::endl;
            // Критерий останова
            if (L < U_RELATION_CONST) {
                break;
            }
            x = shiftBaseNodesGrid(x, localExtremums);
            ++it;
        }
        print::vector(x, true);
    }

    // Нужно ли пересчитывать сетку (1 и каждую 10 итерацию)
    bool mustRecalculateGrid(size_t it) {
        return it % 10 == 0;
    }

    // Сформировать сетку дополнительных точек профиля погрешности на основе сдвига узлов
    std::vector<NodeInfo> formAdditionalGrid(const std::vector<NodeInfo>& firstIterationGrid, const BmpVector& shiftedBaseGrid) {
        std::vector<NodeInfo> result{ firstIterationGrid.front() };
        size_t count = 1;
        for (size_t i = 1; i < firstIterationGrid.size() - 1; ++i) {
            if (firstIterationGrid.at(i).value > shiftedBaseGrid.at(count)) {
                result.emplace_back(NodeInfo{ shiftedBaseGrid.at(count), true });
                ++count;
            }
            result.emplace_back( NodeInfo{ firstIterationGrid.at(i).value, false });
        }
        result.emplace_back(firstIterationGrid.back());
        return result;
    }

    // Сформировать сетку значений функции в дополнительных точках погрешности на основе сдвига узлов
    std::vector<NodeInfo> formPrecisionValues(const BmpVector& x, const std::vector<NodeInfo>& U_prec) {
        std::vector<NodeInfo> result{ U_prec.front() };
        size_t count = 1;
        BmpVector shiftedPrecFuncVales = getPrecisionValue(x);
        for (size_t i = 1; i < U_prec.size() - 1; ++i) {
            // Значение функции в новом базовом узле
            if (U_prec.at(i).value > shiftedPrecFuncVales.at(count)) {
                result.emplace_back(NodeInfo{ shiftedPrecFuncVales.at(count), true});
                ++count;
            }
            result.emplace_back(NodeInfo{ U_prec.at(i).value, false });
        }
        result.emplace_back(U_prec.back());
        return result;
    }

    BmpVector getInversedCoord(const BmpVector& x, const BmpReal P) {
        BmpVector result;
        for (const auto& it : x) {
            result.emplace_back(BmpReal(1) / pow(it, P));
        }
        return result;
    }

    // Тест на функции ФД
    void testFDFunction() {
        const BmpReal N = 1; // 2+2 коэффициента 
        //const BmpReal N = 3; // 4+4 коэффициента 
        BmpVector x = getLinTrigDistribution(N);
        x.emplace_back(BmpReal(1) / 1e-7); // т. бесконечность
        // Значение константы отношения u_max / u_min - 1
        const BmpReal U_RELATION_CONST = 0.1; // Для ФФД 
        size_t it = 1;
        // TODO: Сделать без пересчета значений функции на каждой итерации
        std::vector<NodeInfo> additionalDots = getAdditionalDotsWithInfo(x);
        std::vector<NodeInfo> U_prec = getPrecisionValueWithInfo(additionalDots);
        const size_t MAX_ITERATION_COUNT = 50;
        BmpReal L_prev = 1e10;
        while (it < MAX_ITERATION_COUNT) {
            std::cout << std::endl << std::endl;
            std::cout << "Current iteration:" << it << std::endl;
            std::cout << "Current grid x: " << std::endl;
            print::vector(x); std::cout << std::endl << std::endl;
            if (mustRecalculateGrid(it)) {
                additionalDots = getAdditionalDotsWithInfo(x);
                U_prec = getPrecisionValueWithInfo(additionalDots);
            } else if (it > 1) {
                additionalDots = formAdditionalGrid(additionalDots, x);
                U_prec = formPrecisionValues(x, U_prec);
            }
            /*
            BmpVector additionalDots = getAdditionalDots(x);
            BmpVector U_prec = getPrecisionValue(additionalDots);
            */
            BmpVector U = approximateRationalFd(x, convNodeInfoToVec(additionalDots), convNodeInfoToVec(U_prec));
            std::vector<NodeInfo> delta = getErrorProfile(U_prec, U);
            // Получаем экстремумы
            // BmpVector localExtremums = localInterpolationExtremums(delta, 2 * N + 1);
            BmpVector localExtremums = localInterpolationExtremumsWithIndependentAddDots(delta, 2 *(N + 1));
            std::cout << "Current extremums: " << std::endl;
            print::vector(localExtremums); std::cout << std::endl;
            BmpVector moduleLocalExtremums = abs(localExtremums);
            // Диагностика
            BmpReal L = *std::max_element(moduleLocalExtremums.begin(), moduleLocalExtremums.end()) / *std::min_element(moduleLocalExtremums.begin(), moduleLocalExtremums.end()) - 1;
            std::cout << "Current L: " << L << std::endl;
            // Критерий останова
            // Если достигли заданной точности
            // Если диагностика на текущей итерации превышает значения на предыдущей, выходим 
            if ((L < U_RELATION_CONST) || (it > 1 && L > L_prev)) {
                break;
            }
            L_prev = L;
            // Сдвиг узлов производится для координат 1/(x^2)
            BmpVector xInversed = getInversedCoord(x, P);
            std::reverse(xInversed.begin(), xInversed.end());
            std::reverse(localExtremums.begin(), localExtremums.end());
            xInversed = shiftBaseNodesGrid(xInversed, localExtremums);
            std::reverse(xInversed.begin(), xInversed.end());
            x = getInversedCoord(xInversed, 1.0 / P);
            ++it;
        }
        print::vector(x, true);
    }

}

TEST_CASE("check") {
   //testZeroFunction();
   testFDFunction();
}