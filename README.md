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
    * AsymptoticSeries.cpp - вычисление асимптотическим рядом
    * FullyConvergedSeries.cpp - вычисление всюду сходящимся рядом
    * ZetaFunction.cpp - значения дзета-функции Римана

* fd-half
    * example-asym-series-calculate.cpp - вычисление ФД асимптотическим рядом на границе x_min
    * example-fcs-calculate.cpp - вычисление ФД всюду сходящимся рядом в точке x = 0
    * example-half-integer.cpp - вычисление ФД квадратурами в промежутке 0 <= x <= x_min