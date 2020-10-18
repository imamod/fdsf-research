#include "Common.h"
#include "Grid.h"
#include "MatrixUtils.h"
#include "Fdsf.h"
#include "FdIndex.h"
#include "Gamma.h"
#include <algorithm>

/* Множители имеют одинаковый знак */
class MultipliersHaveTheSameSign : public std::exception {
    public:
        const char* what() const {
            return "Multipliers have the same sign";
        }
};

namespace {

    // Возвращает вектор абсолютных значений
    BmpVector abs(const BmpVector& v) {
        BmpVector absoluteValue;
        for (auto it : v) {
            absoluteValue.push_back(std::abs(it));
        }
        return absoluteValue;
    }

    // Отладка автоматического выравниваия погрешности
    const BmpReal one = BmpReal(1);
    const BmpReal num2 = BmpReal(2);
    const BmpReal x_star = BmpReal(4);
    // 
    BmpReal m_alphaRight = 2 / (2 + pi());
    BmpReal P(2);

    // Установка веса линейно-тригонометрической сетки и начального приближения в зависимости от k
    void setModuleVariables(BmpReal k) {
        if (fdsf::index::M3_HALF == k || fdsf::index::P7_HALF == k) {
            m_alphaRight = 0.6;
            P = 1.5;
        } else {
            m_alphaRight = 2 / (2 + pi());
            P = 2.0;
        }
    }

    // Сетка узлов для тестирования двухкусковой формулы
    class GridIhalf : public Grid {
        public:
            GridIhalf(size_t N_base, size_t addNCount)
                : Grid(N_base, addNCount) {}

            // Установить сетку базовых узлов по линейно-тригонометрическому закону
            virtual void setLinearTrigonometricGrid() {
                const BmpReal alpha = 2 / (2 + pi());
                const BmpReal y_star = BmpReal(log(1 + exp(x_star)));
                BmpReal baseSize = BmpReal(2 * baseNCount() - 1);
                // Задаются базовые узлы интерполяции
                BmpVector baseGrid;
                for (size_t j = 1; j <= baseSize; j++) {
                    baseGrid.push_back(y_star / num2 * (num2 * alpha * j / baseSize
                        + (one - alpha) * (one - cos(pi() * j / baseSize))));
                }
                setBaseGrid(baseGrid);
            }

            // Задает правую линейно-тригонометрическую сетку базовых узлов справа
            virtual void setLinearTrigonometricGridRight() {
                const BmpReal y_star = x_star;
                BmpReal baseSize = BmpReal(2 * baseNCount());

                const BmpReal y_star_new = BmpReal(1) / pow(y_star, P);
                // Задаются базовые узлы интерполяции
                BmpVector baseGrid;
                // Задаются базовые узлы интерполяции
                for (size_t j = 1; j <= baseSize; j++) {
                    /* Усиленная линейная часть */
                    BmpReal linearPart = num2 * (one - m_alphaRight) * j / baseSize;
                    BmpReal trigonometricPart = m_alphaRight * (one - cos(pi() * j / baseSize));
                    baseGrid.push_back(pow(y_star_new / num2 * (linearPart + trigonometricPart), 1 / P));
                }

                // Разворачиваем y
                std::reverse(baseGrid.begin(), baseGrid.end());
                changeGrid(baseGrid);
                setBaseGrid(baseGrid);
            }
    };

    // Вспомогательная структура с информацией о узле интерполяции
    struct NodeInfo {
        BmpReal value;    // значение точки
        bool isBaseNode;  // является ли базовым
    };

    void printBaseNodes(const std::vector<NodeInfo>& nodes) {
        print::string("base node: ");
        for (auto it : nodes) {
            if (it.isBaseNode) {
                std::cout << "    " << it.value << std::endl;
            }
        }
    }

    BmpVector getBaseNodesError(const std::vector<NodeInfo>& nodes) {
        BmpVector baseNodes;
        for (auto it : nodes) {
            if (it.isBaseNode) {
                baseNodes.emplace_back(it.value);
            }
        }
        return baseNodes;
    }

    void printNodes(const std::vector<NodeInfo>& nodes) {
        print::string("nodes: ");
        for (auto it : nodes) {
            std::cout << "    " << it.value << std::endl;
        }
    }

    BmpVector convNodeInfoToVec(const std::vector<NodeInfo>& v) {
        BmpVector result;
        for (auto it : v) {
            result.emplace_back(it.value);
        }
        return result;
    }

    const size_t ADDITIONAL_DOTS_COUNT = 20;

    /* Алгоритм выравнивания экстремумов */

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
                throw MultipliersHaveTheSameSign();
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
        const BmpReal A = 0.2; // Release
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

    /****** Методический расчет*******/

    // Получает локальный экстремум на заданном промежутке
    template <typename It>
    BmpReal getLocalExtremum(It first, It last, bool isPositive) {
        return isPositive ? *std::max_element(first, last) : *std::min_element(first, last);
    }

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

    /**
     * Аппроксимация полиномом
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
    BmpVector getUniformDistributionDots(size_t N) {
        BmpVector x;
        for (size_t n = 0; n <= N; ++n) {
            x.push_back(BmpReal(2 * n - N) / N);
        }
        return x;
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

    void testZeroFunction() {
        const size_t N = 3;
        BmpVector x = getUniformDistributionDots(N);
        // Значение константы отношения u_max / u_min - 1
        const BmpReal U_RELATION_CONST = 0.01; // Для модельных примеров
        size_t it = 1;
        const size_t MAX_ITERATION_COUNT = 50;
        while (it < MAX_ITERATION_COUNT) {
            std::cout << std::endl << std::endl;
            std::cout << "Current iteration:" << it << std::endl;
            print::string("Current grid x: ");
            print::vector(x); std::cout << std::endl << std::endl;
            BmpVector additionalDots = getAdditionalDots(x);
            // Аппроксимирующий многочлен (для U = 0 погрешности - сами значения аппроксимации) 
            BmpVector U = approximatePolinom(x, additionalDots);
            BmpVector delta = U;
            // Получаем экстремумы
            BmpVector localExtremums = localInterpolationExtremums(delta, N);
            print::string("Current extremums: ");
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

    /************Правая аппроксимация ФД*************/

    // Литр справа, Т - число экстремумов
    BmpVector getLinTrigDistributionRight(size_t N) {
        GridIhalf grid(N + 1, ADDITIONAL_DOTS_COUNT);
        grid.setLinearTrigonometricGridRight();
        return grid.base();
    }

    std::function<BmpReal(BmpReal)> fdFunc(BmpReal k) {
        std::function<BmpReal(BmpReal)> f;
        if (fdsf::index::M3_HALF == k) {
            f = fdsf::fd_m3half;
        } else if (fdsf::index::M1_HALF == k) {
            f = fdsf::fd_m1half;
        } else if (fdsf::index::P1_HALF == k) {
            f = fdsf::fd_1half;
        } else if (fdsf::index::P3_HALF == k) {
            f = fdsf::fd_3half;
        } else if (fdsf::index::P5_HALF == k) {
            f = fdsf::fd_5half;
        } else if (fdsf::index::P7_HALF == k) {
            f = fdsf::fd_7half;
        }
        return f;
    }

    BmpVector getPrecisionValue(const BmpVector& x, BmpReal k) {
        BmpVector result;
        auto f = fdFunc(k);
        for (auto it : x) {
            result.emplace_back(f(it));
        }
        return result;
    }

    std::vector<NodeInfo> getPrecisionValueWithInfo(BmpReal k, const std::vector<NodeInfo>& x) {
        std::vector<NodeInfo> result;
        auto f = fdFunc(k);
        for (auto it : x) {
            result.emplace_back(NodeInfo{ f(it.value), it.isBaseNode });
        }
        print::string("result size " + std::to_string(result.size()));
        return result;
    }


    //// Аппроксимирующая формула 1 
    BmpVector getApproximationFormByFDValue(BmpReal k, const BmpVector& y, const BmpVector& v) {
        const BmpReal C1 = (k + 1) * k * pow(pi(), 2) / 6;
        BmpVector result;
        for (size_t i = 0; i < v.size(); ++i) {
            result.emplace_back(((v[i] * (k + 1) / pow(y[i], k + 1) - 1) * pow(y[i], 2)) / C1);
        }
        return result;
    }

    BmpVector getFDValueFromApproximation(BmpReal k, const BmpVector& Y, const BmpVector& v) {
        BmpVector U;
        const BmpReal C1 = (k + 1) * k * pow(pi(), 2) / 6;
        for (size_t i = 0; i < v.size(); ++i) {
            U.emplace_back((v.at(i) * C1 / pow(Y.at(i), 2) + 1) * (pow(Y.at(i), k + 1) / (k + 1)));
        }
        return U;
    }

    // Аппроксимирующая форма 2
    BmpVector getApproximationFormByFDValue2(BmpReal k, const BmpVector& y, const BmpVector& v) {
        const BmpReal C2 = (k + 1) * pow(pi(), 2) / 3;
        BmpVector result;
        for (size_t i = 0; i < v.size(); ++i) {
            result.emplace_back(((pow(v[i] * (k + 1) / y[i], BmpReal(2.0) / k) - pow(y[i], 2))) / C2);
        }
        return result;
    }

    BmpVector getFDValueFromApproximation2(BmpReal k, const BmpVector& Y, const BmpVector& v) {
        BmpVector U;
        const BmpReal C2 = (k + 1) * pow(pi(), 2) / 3;
        for (size_t i = 0; i < v.size(); ++i) {
            U.emplace_back(pow(v.at(i) * C2 + pow(Y.at(i), 2), k / BmpReal(2.0)) * Y.at(i) / (k + 1));
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

    BmpVector localInterpolationExtremumsWithIndependentAddDots(const std::vector<NodeInfo>& errorProfile, size_t N) {
        BmpVector localExtremums;
        auto isBaseNode = [](const NodeInfo& x) {
            return x.isBaseNode;
        };
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

    BmpVector approximateRationalFd(BmpReal k, const BmpVector& x, const BmpVector& X, const BmpVector& U_additional) {

        BmpVector U_base = getPrecisionValue(x, k);
        const size_t baseSize = x.size()-1; // т. БЕсконечность в аппроксимации не участвует
        const size_t N = baseSize / 2;
        BmpVector y0 = getI0Values(x);
        BmpVector Y = getI0Values(X);
        //// Первая форма аппроксимации
        //BmpVector z = getApproximationFormByFDValue(k, y0, U_base);
        //BmpVector I = getApproximationFormByFDValue(k, Y, U_additional);
        //// Вторая форма аппроксимации
        BmpVector z = getApproximationFormByFDValue2(k, y0, U_base);
        BmpVector I = getApproximationFormByFDValue2(k, Y, U_additional);

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

        print::string("Аппроксимирующие коэффициенты в числителе:");
        print::vector(a);
        print::string("Аппроксимирующие коэффициенты в знаменателе:");
        print::vector(b);

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

        //// Первая форма аппроксимации
        //return getFDValueFromApproximation(k, Y, F_additional);
        //// Вторая форма аппроксимации
        return getFDValueFromApproximation2(k, Y, F_additional);
    }

    // Для литр, доп точки нужно получать точки по обратной величине
    std::vector<NodeInfo> getAdditionalDotsWithInfo(const BmpVector& x) {
        BmpVector inverseValue;
        for (const auto& it : x) {
            inverseValue.emplace_back(BmpReal(1) / pow( it, P));
        }
        print::string("x size " + std::to_string(x.size()));
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
        return additionalDots;
    }

    std::vector<NodeInfo> getErrorProfile(const std::vector<NodeInfo>& U_prec, const BmpVector& U) {
        std::vector<NodeInfo> delta;
        for (size_t i = 0; i < U_prec.size(); ++i) {
            delta.emplace_back(NodeInfo{ U[i] / U_prec[i].value - 1, U_prec[i].isBaseNode });
        }
        return delta;
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
    std::vector<NodeInfo> formPrecisionValues(BmpReal k, const BmpVector& x, const std::vector<NodeInfo>& U_prec) {
        std::vector<NodeInfo> result{ U_prec.front() };
        size_t count = 1;
        BmpVector shiftedPrecFuncVales = getPrecisionValue(x, k);
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

    // Получить значение диагностики L
    BmpReal getLDiagnostic(const BmpVector& extremums) {
        return *std::max_element(extremums.begin(), extremums.end()) / *std::min_element(extremums.begin(), extremums.end()) - 1;
    }

    BmpReal getSKO(const BmpVector& x) {
        BmpReal result = 0;
        for (auto it : x) {
            result += it * it;
        }
        return sqrt(result / x.size());
    }

    // Получение СКО-диагностики профиля и узлов
    BmpReal getSKODiagnostic(const BmpVector& base, const BmpVector& extr) {
        BmpReal nom = getSKO(base);
        BmpReal denom = getSKO(extr);
        print::string("base SKO: " + std::to_string(nom));
        print::string("extr SKO: " + std::to_string(denom));
        return nom / denom;
    }

    // Тест на функции ФД
    void testFDFunction(BmpReal k) {
        setModuleVariables(k);
       // const size_t N = 0; // 1+1коэффициент
       // const size_t N = 1; // 2+2 коэффициента
       // const size_t N = 2; // 3+3 коэффициента
        const size_t N = 3; // 4+4 коэффициента Release
        BmpVector x = getLinTrigDistributionRight(N);
        x.emplace_back(BmpReal(1) / 1e-7); // т. бесконечность
        // Значение константы отношения u_max / u_min - 1
        const BmpReal U_RELATION_CONST = 0.1; // Для ФФД 
        size_t it = 1;
        std::vector<NodeInfo> additionalDots = getAdditionalDotsWithInfo(x);
        std::vector<NodeInfo> U_prec = getPrecisionValueWithInfo(k, additionalDots);
        const size_t MAX_ITERATION_COUNT = 50;
        BmpReal L_prev = 1e10;
        while (it < MAX_ITERATION_COUNT) {
            std::cout << std::endl << std::endl;
            print::string("Current iteration:" + std::to_string(it));
            print::string("Current grid x: ");
            print::vector(x); std::cout << std::endl << std::endl;
            if (mustRecalculateGrid(it)) {
                additionalDots = getAdditionalDotsWithInfo(x);
                U_prec = getPrecisionValueWithInfo(k, additionalDots);
            } else if (it > 1) {
                additionalDots = formAdditionalGrid(additionalDots, x);
                U_prec = formPrecisionValues(k, x, U_prec);
            }
            /*
            BmpVector additionalDots = getAdditionalDots(x);
            BmpVector U_prec = getPrecisionValue(additionalDots);
            */
            BmpVector U = approximateRationalFd(k, x, convNodeInfoToVec(additionalDots), convNodeInfoToVec(U_prec));
            std::vector<NodeInfo> delta = getErrorProfile(U_prec, U);
            printBaseNodes(delta);
            // Получаем экстремумы
            // BmpVector localExtremums = localInterpolationExtremums(delta, 2 * N + 1);
            BmpVector localExtremums = localInterpolationExtremumsWithIndependentAddDots(delta, 2 *(N + 1));
            print::string("Current extremums: ");
            print::vector(localExtremums); std::cout << std::endl;
            BmpVector moduleLocalExtremums = abs(localExtremums);
            // Диагностика СКО
            BmpReal SKO = getSKODiagnostic(getBaseNodesError(delta), localExtremums);
            print::string("Current SKO: " + std::to_string(SKO));
            // Диагностика L
            BmpReal L = getLDiagnostic(moduleLocalExtremums);
            std::cout << "Current L: " << L << std::endl;
            // Критерий останова
            // Если достигли заданной точности
            // Если диагностика на текущей итерации превышает значения на предыдущей, выходим 
            if ((L < U_RELATION_CONST)) { // (it > 1 && L > L_prev)
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


    /************Левая аппроксимация ФД*************/

    // Возвращает линейно-тригонометрическое распределение узлов для левой аппроксимации
    BmpVector getLinTrigDistributionLeft(size_t N) {
        GridIhalf grid(N + 1, ADDITIONAL_DOTS_COUNT);
        grid.setLinearTrigonometricGrid();
        return grid.base();
    }

    // Для литр, доп точки нужно получать точки по обратной величине
    std::vector<NodeInfo> getLeftAdditionalDotsWithInfo(const BmpVector& _x) {
        std::cout << "x size " << _x.size() << std::endl;
        // Для левой асимптотики расчет ведется на промежутке (0;x*],
        // поэтому добавляем фиктивную границу 0
        BmpVector x{ 0 };
        x.insert(x.end(), _x.begin(), _x.end());
        std::vector<NodeInfo> additionalDots;
        for (size_t i = 0; i < x.size() - 1; ++i) {
            BmpReal delta = (x.at(i + 1) - x.at(i)) / ADDITIONAL_DOTS_COUNT;
            for (size_t j = 1; j < ADDITIONAL_DOTS_COUNT; ++j) {
                additionalDots.emplace_back(NodeInfo{ x.at(i) + j * delta, false });
            }
            additionalDots.emplace_back(NodeInfo{ x.at(i + 1), true });
        }
        std::cout << "additionalDots size " << additionalDots.size() << std::endl;
        return additionalDots;
    }


    //// Левая аппроксимирующая формула
    BmpVector getApproximationFormByFDValueLeft(BmpReal k, const BmpVector& y, const BmpVector& v) {
        const BmpReal gammaValue = factorial(k);
        BmpVector result;
        for (size_t i = 0; i < v.size(); ++i) {
            result.emplace_back(pow(v.at(i) / (gammaValue * y.at(i)), BmpReal(1.0)/k));
        }
        return result;
    }

    BmpVector getFDValueFromApproximationLeft(BmpReal k, const BmpVector& Y, const BmpVector& v) {
        BmpVector U;
        const BmpReal gammaValue = factorial(k);
        for (size_t i = 0; i < v.size(); ++i) {
            U.emplace_back(pow(v.at(i), k) * Y.at(i) * gammaValue);
        }
        return U;
    }

    // Левая аппроксимирующая формула TODO:
    BmpVector approximateRationalFdLeft(BmpReal k, const BmpVector& x, const BmpVector& X, const BmpVector& U_additional) {

        BmpVector U_base = getPrecisionValue(x, k);
        const size_t baseSize = x.size();
        const size_t N = (baseSize-1) / 2;
        BmpVector y0 = getI0Values(x);
        BmpVector Y = getI0Values(X);
        //// Первая форма аппроксимации
        BmpVector z = getApproximationFormByFDValueLeft(k, y0, U_base);
        BmpVector I = getApproximationFormByFDValueLeft(k, Y, U_additional);

        // Правая часть
        BmpVector rightPart;
        for (auto it : z) {
            rightPart.emplace_back(it - 1);
        }
        // Левая часть
        BmpMatrix A(baseSize, BmpVector(baseSize, 0));
        for (size_t i = 0; i < baseSize; ++i) {
            for (size_t j = 0; j < baseSize; ++j) {
                if (j < N + 1) {
                    A.at(i).at(j) = pow(y0.at(i), BmpReal(j + 1));
                } else {
                    A.at(i).at(j) = -z.at(i) * pow(y0.at(i), BmpReal(j - N));
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
            if (j < N + 1) {
                a.emplace_back(E.at(j));
            } else if (j <= baseSize) {
                b.emplace_back(E.at(j));
            }
        }

        print::string("Аппроксимирующие коэффициенты в числителе:");
        print::vector(a);
        print::string("Аппроксимирующие коэффициенты в знаменателе:");
        print::vector(b);

        // Вычисление аппроксимирующей величины
        BmpVector F_additional;
        for (auto y : Y) {
            BmpReal nom = 0, denom = 0;
            for (size_t n = 1; n <= N + 1; ++n) {
                nom += a.at(n - 1) * pow(y, BmpReal(n));
            }
            for (size_t m = 1; m <= N; ++m) {
                denom += b.at(m - 1) * pow(y, BmpReal(m));
            }
            F_additional.emplace_back(BmpReal(1.0 + nom) / (1.0 + denom));
        }
        // Восстановления значения функции по аппроксимации
        return getFDValueFromApproximationLeft(k, Y, F_additional);
    }

    void testFDLeftApproximation(BmpReal k) {
       // const size_t N = 1; // 2+1 коэффициент
       // const size_t N = 2; // 3+2 коэффициента
        const size_t N = 3; // 4+3 коэффициента Release
       // const size_t N = 4; // 5+4 коэффициента
       // const size_t N = 5; // 6+5 коэффициента
       // const size_t N = 6; // 7+6 коэффициента
        BmpVector x = getLinTrigDistributionLeft(N);
        // Значение константы отношения u_max / u_min - 1
        const BmpReal U_RELATION_CONST = 0.1; // Для ФФД 
        size_t it = 1;
        std::vector<NodeInfo> additionalDots = getLeftAdditionalDotsWithInfo(x);
        std::vector<NodeInfo> U_prec = getPrecisionValueWithInfo(k, additionalDots);
        const size_t MAX_ITERATION_COUNT = 50;
        BmpReal L_prev = 1e10;
        while (it < MAX_ITERATION_COUNT) {
            std::cout << std::endl << std::endl;
            print::string("Current iteration:" + std::to_string(it));
            print::string("Current grid x: ");
            print::vector(x); std::cout << std::endl << std::endl;
            if (mustRecalculateGrid(it)) {
                additionalDots = getLeftAdditionalDotsWithInfo(x);
                U_prec = getPrecisionValueWithInfo(k, additionalDots);
            } else if (it > 1) {
                // Добавляем нулевую, фиктивную точку
                x.insert(x.begin(), 0);
                additionalDots = formAdditionalGrid(additionalDots, x);
                U_prec = formPrecisionValues(k, x, U_prec);
                // Удаляем нулевую, фиктивную точку
                x.erase(x.begin());
            }
            BmpVector U = approximateRationalFdLeft(k, x, convNodeInfoToVec(additionalDots), convNodeInfoToVec(U_prec));
            std::vector<NodeInfo> delta = getErrorProfile(U_prec, U);
            printBaseNodes(delta);
            printNodes(delta);
            // Получаем экстремумы
            BmpVector localExtremums = localInterpolationExtremumsWithIndependentAddDots(delta, 2 * N + 1);
            print::string("Current extremums: ");
            print::vector(localExtremums); std::cout << std::endl;
            BmpVector moduleLocalExtremums = abs(localExtremums);
            // Диагностика L
            BmpReal L = getLDiagnostic(moduleLocalExtremums);
            std::cout << "Current L: " << L << std::endl;
            // Критерий останова
            // Если достигли заданной точности
            // Если диагностика на текущей итерации превышает значения на предыдущей, выходим 
            if ((L < U_RELATION_CONST)) { //|| (it > 1) && L > L_prev)) { // Временно отключили диагностику L
                break;
            }
            L_prev = L;
            // Добавляем нулевую, фиктивную точку
            x.insert(x.begin(), 0);
            // Сдвиг узлов производится для координат y
            x = shiftBaseNodesGrid(x, localExtremums);
            // Удаляем нулевую, фиктивную точку
            x.erase(x.begin());
            ++it;
        }
        print::vector(x, true);
    }
}

TEST_CASE("zero") {
    testZeroFunction();
}

TEST_CASE("fd_right") {
    SECTION("km32") {
        testFDFunction(fdsf::index::M3_HALF);
    }
    SECTION("km12") {
        testFDFunction(fdsf::index::M1_HALF);
    }
    SECTION("k12") {
        testFDFunction(fdsf::index::P1_HALF);
    }
    SECTION("k32") {
        testFDFunction(fdsf::index::P3_HALF);
    }
    SECTION("k52") {
        testFDFunction(fdsf::index::P5_HALF);
    }
    SECTION("k72") {
        testFDFunction(fdsf::index::P7_HALF);
    }
}

TEST_CASE("fd_left") {
    SECTION("km32") {
        testFDLeftApproximation(fdsf::index::M3_HALF);
    }
    SECTION("km12") {
        testFDLeftApproximation(fdsf::index::M1_HALF);
    }
    SECTION("k12") {
        testFDLeftApproximation(fdsf::index::P1_HALF);
    }
    SECTION("k32") {
        testFDLeftApproximation(fdsf::index::P3_HALF);
    }
    SECTION("k52") {
        testFDLeftApproximation(fdsf::index::P5_HALF);
    }
    SECTION("k72") {
        testFDLeftApproximation(fdsf::index::P7_HALF);
    }
}
