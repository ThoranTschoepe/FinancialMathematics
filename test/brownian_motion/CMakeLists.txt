set(EXAMPLE_NAME test-${NAME}-brownian_motion)

project(${EXAMPLE_NAME})

add_executable(${EXAMPLE_NAME} main.cpp)

target_link_libraries(${EXAMPLE_NAME} ${NAME} gtest gtest_main)