#pragma once

#include <exception>

/* ������� �� ����������� */
class FunctionNotImplemented : public std::exception {
    public:
        virtual const char* what() const {
            return "Function not implemented";
        }
};

/* ���������� ������ ���������� ������� */
class UnknownCalculationMethod : public std::exception {
    public:
        virtual const char* what() const {
            return "Unknown calculation method of FD-function";
        }
};

/* ���������������� ������� �� */
class UnsupportedFdFunction : public std::exception {
    public:
        virtual const char* what() const {
            return "Unsupported FD-function index";
        }
};
