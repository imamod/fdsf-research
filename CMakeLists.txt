cmake_minimum_required(VERSION 2.8)

project(fdsf-research)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

include_directories("src/include")

file(GLOB LIB_SOURCES "src/*.cpp")

set (SOURCE_FILES ${LIB_SOURCES})
add_library (fdsf-research STATIC ${SOURCE_FILES})
