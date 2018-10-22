cmake_minimum_required(VERSION 2.8)
project(example-exponent-convergence)

set(CATCH_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../3dparty)
message(STATUS "CATCH_INCLUDE_DIR: " ${CATCH_INCLUDE_DIR})

set(SRC example-exponent-convergence.cpp)
add_executable(${SRC} ${project})

target_include_directories(example-exponent-convergence PRIVATE ${CATCH_INCLUDE_DIR})

if (MSVC)
    target_link_libraries(${testcase} fdsf.lib Utils.lib)
else()
    target_link_libraries(${testcase} Utils.a fdsf.a stdc++fs)
endif()
