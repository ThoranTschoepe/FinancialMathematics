if (PROJECT_IS_TOP_LEVEL)
    include(FetchContent)

    FetchContent_Declare(
            Eigen
            GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
            GIT_TAG master
            GIT_SHALLOW TRUE
            GIT_PROGRESS TRUE)

    set(EIGEN_BUILD_DOC OFF)
    set(BUILD_TESTING OFF)
    set(EIGEN_BUILD_PKGCONFIG OFF)
    set( OFF)

    FetchContent_Declare(
            googletest
            URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.tar.gz
            DOWNLOAD_EXTRACT_TIMESTAMP 1
    )
    FetchContent_Declare(
            googlebenchmark
            URL https://github.com/google/benchmark/archive/refs/tags/v1.8.3.tar.gz
            DOWNLOAD_EXTRACT_TIMESTAMP 1
    )

    set (BENCHMARK_ENABLE_INSTALL OFF)
    set (BENCHMARK_ENABLE_TESTING OFF)
    FetchContent_MakeAvailable(Eigen googletest googlebenchmark)
endif ()