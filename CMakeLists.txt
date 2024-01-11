cmake_minimum_required(VERSION 3.27)

#set name
set(NAME FinancialMathematics)

project(${NAME})

set(CMAKE_CXX_STANDARD 20)

set(BUILD_EXAMPLES_${NAME} ON CACHE BOOL "BUILD_EXAMPLES" FORCE)
set(BUILD_TESTS_${NAME} ON CACHE BOOL "BUILD_TESTS" FORCE)
set(BUILD_BENCHMARK_${NAME} ON CACHE BOOL "BUILD_BENCHMARK" FORCE)

include(FetchDependencies.cmake)

add_library(${NAME}
        src/brownian_motion.cpp
        src/random_utils.cpp
        src/mc.cpp
        src/AssetPaths.cpp
        src/black_scholes.cpp
)

target_include_directories(${NAME}
        PUBLIC include
)

find_package(OpenMP REQUIRED)

set(EIGEN_USE_BLAS ON)

target_link_libraries(${NAME} Eigen3::Eigen OpenMP::OpenMP_CXX)

if (PROJECT_IS_TOP_LEVEL)
    if (BUILD_EXAMPLES_${NAME})
        add_subdirectory(examples)
    endif ()
    if (BUILD_BENCHMARK_${NAME})
        add_subdirectory(benchmark)
    endif ()
    if (BUILD_TESTS_${NAME})
        add_subdirectory(test)
    endif ()
endif ()