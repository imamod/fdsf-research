#pragma once

#include <exception>

/* Функция не реализована */
class FunctionNotImplemented : public std::exception {
    public:
        virtual const char* what() const {
            return "Function not implemented";
        }
};

/* Неизестный способ вычислений функции */
class UnknownCalculationMethod : public std::exception {
    public:
        virtual const char* what() const {
            return "Unknown calculation method of FD-function";
        }
};

/* Неподдерживаемая функция ФД */
class UnsupportedFdFunction : public std::exception {
    public:
        virtual const char* what() const {
            return "Unsupported FD-function index";
        }
};
