#pragma once
#define LOGGER_SUPPORTED
#include <string>
#include <fstream>

const std::string LOG_NAME = "Fdsf.log";

class Logger {
    public:
        explicit Logger(const std::string& text)
        : m_text(text)
        , m_file(LOG_NAME, std::ios::app){
#ifdef LOGGER_SUPPORTED
            m_file << "[BEGIN].." << text << std::endl;
#endif
        }

        ~Logger() {
#ifdef LOGGER_SUPPORTED
            m_file << "[END]...." << m_text << std::endl;
            m_file.close();
#endif
        }

        // �������� � ��� ����������
        void info(const std::string& text) {
            write("INFO", text);
        }

        // �������� � ��� �������� �� ������
        void error(const std::string& text) {
            write("ERROR", text);
        }

    private:
        // �����
        std::string m_text;
        // ����
        std::fstream m_file;

        void write(const std::string& tag, const std::string& text) {
#ifdef LOGGER_SUPPORTED
            m_file << "[" << tag << "]...." << text << std::endl;
#endif
        }
};
