cmake_minimum_required(VERSION 2.8)
#############################################################################
# Catch library with the main function to speed up build
#############################################################################
project(Test)

# https://stackoverflow.com/questions/2368811/how-to-set-warning-level-in-cmake
if(MSVC)
    # Force to always compile with W4
    if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
        string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
    endif()

    # Disable warning C4389: '==': signed/unsigned mismatch
    # Disable warning C4309: 'static_cast': truncation of constant value
    # Disable warning C4566: character represented by universal-character-name '\uFF01' cannot be represented in the current code page (1252)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4389 /wd4309 /wd4566")
endif()

# Prepare "Catch" library for other executables
# https://github.com/catchorg/Catch2/blob/master/docs/build-systems.md
set(CATCH_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../3dparty)
message(STATUS "CATCH_INCLUDE_DIR: " ${CATCH_INCLUDE_DIR})
#add_library(Catch INTERFACE)
#target_include_directories(Catch INTERFACE ${CATCH_INCLUDE_DIR})

#############################################################################
# Static library Utils.lib
#############################################################################
set(UTILS_FILES src/Common.cpp src/Common.h)
add_library(Utils STATIC ${UTILS_FILES})
if (NOT MSVC)
    set_target_properties(Utils PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wno-deprecated -Wno-float-equal")
endif()
target_include_directories(Utils PRIVATE ${CATCH_INCLUDE_DIR})
target_include_directories(Utils PUBLIC "../../../../boost_1_60_0")
target_include_directories(Utils PRIVATE "../src/include")
target_include_directories(Utils PRIVATE "../include")
target_include_directories(Utils PUBLIC "../../Utils/src")
target_include_directories(Utils PRIVATE "src")
if ("${CMAKE_BUILD_TYPE}" STREQUAL "DebugMP")
   add_definitions(-DHIGH_PRECISION)
endif ()


#############################################################################
# one executable for each unit test file
#############################################################################

set(FDSF_LIB ${CMAKE_CURRENT_SOURCE_DIR}/../build/${CMAKE_GENERATOR_PLATFORM})
message(STATUS "FDSF_LIB: " ${FDSF_LIB})
if (MCVC)
    set(UTILS_LIB ${CMAKE_CURRENT_SOURCE_DIR}/build/Debug)
else ()
    set(UTILS_LIB ${CMAKE_CURRENT_SOURCE_DIR}/build)
endif()
message(STATUS "UTILS_LIB: " ${UTILS_LIB})

file(GLOB projects "src/example-*.cpp"
                   "src/test-*.cpp"
                   "src/unit-*.cpp"
                   "src/util-*.cpp")
foreach(project ${projects})
    get_filename_component(file_basename ${project} NAME_WE)
    string(REGEX REPLACE "unit-([^$]+)" "unit-\\1" testcase ${file_basename})
    link_directories(${FDSF_LIB} ${UTILS_LIB})

    add_executable(${testcase} ${project})
    add_dependencies(${testcase} Utils)
    if (MSVC)
        set(CMAKE_CONFIGURATION_TYPES Debug Release DebugMP ReleaseMP)
        set(CMAKE_CONFIGURATION_TYPES "${CMAKE_CONFIGURATION_TYPES}" CACHE STRING
            "Reset the configurations to what we need"
            FORCE)
        cmake_policy(SET CMP0043 NEW)
        set_target_properties(${testcase} PROPERTIES
            CXX_STANDARD 14
            CXX_STANDARD_REQUIRED ON
            COMPILE_DEFINITIONS "$<$<CXX_COMPILER_ID:MSVC>:_SCL_SECURE_NO_WARNINGS>"
            COMPILE_OPTIONS "$<$<CXX_COMPILER_ID:MSVC>:/EHsc;$<$<CONFIG:Release>:/Od>>"
            COMPILE_DEFINITIONS_DEBUGMP "$<$<CXX_COMPILER_ID:MSVC>:HIGH_PRECISION>"
            COMPILE_DEFINITIONS_RELEASEMP "$<$<CXX_COMPILER_ID:MSVC>:HIGH_PRECISION>"
        )
    else()
        set_target_properties(${testcase} PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wno-deprecated -Wno-float-equal")
    endif()

    target_compile_definitions(${testcase} PRIVATE CATCH_CONFIG_FAST_COMPILE)
    target_include_directories(${testcase} PRIVATE ${CATCH_INCLUDE_DIR})
    target_include_directories(${testcase} PRIVATE "src")
    target_include_directories(${testcase} PRIVATE "../src/include")
    target_include_directories(${testcase} PUBLIC "../include")
    target_include_directories(${testcase} PUBLIC "../../Utils/include")
    target_include_directories(${testcase} PUBLIC "../../Utils/src")
    target_include_directories(${testcase} PUBLIC "../../../../boost_1_60_0")
    if (MSVC)
        target_link_libraries(${testcase} fdsf.lib Utils.lib)
    else()
        target_link_libraries(${testcase} Utils.a fdsf.a stdc++fs)
    endif()

endforeach()
