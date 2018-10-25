cmake_minimum_required(VERSION 2.8)

set(project_name Utils)
project(${project_name})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

file(GLOB SOURCES "src/*.cpp")

add_library (${project_name} STATIC ${SOURCES})

target_compile_definitions(${project_name} PRIVATE CATCH_CONFIG_FAST_COMPILE)

set(THIRDPARTY_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../3dparty)
target_include_directories(${project_name} PRIVATE "src/include")
target_include_directories(${project_name} PRIVATE "../fdsf/src/include")
target_include_directories(${project_name} PRIVATE ${THIRDPARTY_INCLUDE_DIR})
