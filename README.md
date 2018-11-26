# fdsf-research

## Модули fdsf-research

* src
   * exp-conv - тесты на экспоненциальную сходимость
   * fdsf - библиотека работы с ФД
   * fd-integer - тесты для ФД целых индексов
   * fd-half - тесты для ФД полуцелых индексов
   * fd-Jmhalf - тесты для интегральной ФД
   * Utils - вспомогательные библиотека
   * high-precision - Wip

## Сборка
### Linux
Запустить скрипт src/build.sh. Соберутся библиотеки Utils, fdsf, тесты для ФД полуцелого индекса и интегральной ФД.

### Windows
Сборка через solution fdsf-research.sln. После сборки solution перейти в директорию с тестами и запустить python-скрипт для требуемых тестов.

## Пояснение модулей

* fdsf
    * AsymptoticSeries.cpp - модуль вычисления асимптотическим рядом
    * FdHalfQuadratures.cpp - модуль вычисления экспоненциально сходящимися квадратурами
    * Fdsf.cpp - модуль вычисления функций ФД целых и полуцелых индексов (-3/2 <= k <= 4)
    * FullyConvergedSeries.cpp - модуль вычисления всюду сходящимся рядом
    * Gamma.cpp - необходимые значения Г-функции
    * ZetaFunction.cpp - значения дзета-функции Римана

* fd-half
    * example-asym-series-calculate.cpp - пример вычисления ФД асимптотическим рядом на границе x_min
    * example-fcs-calculate.cpp - пример вычисления ФД всюду сходящимся рядом в точке x = 0
    * example-half-integer-quadratures.cpp - пример вычисления ФД квадратурами в промежутке 0 <= x <= x_min
    * example-calculate-fd-half.cpp - пример вычисления ФД в произвольном диапазоне

* fd-integer
    * example-fcs-calculate.cpp - пример вычисления ФД всюду сходящимся рядом в точке x = 0
    * example-fd-calculate.cpp - пример вычисления ФД в произвольном диапазоне

* fd-Jmhalf
