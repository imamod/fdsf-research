cmake_minimum_required(VERSION 2.8)

set(project_name fdsf)
project(${project_name})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

include_directories("src/include")

file(GLOB LIB_SOURCES "src/*.cpp")

set (SOURCE_FILES ${LIB_SOURCES})
add_library (${project_name} STATIC ${SOURCE_FILES})
