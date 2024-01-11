set(EXAMPLE_NAME example-${NAME}-vanilla)

project(${EXAMPLE_NAME})

add_executable(${EXAMPLE_NAME} main.cpp)

target_link_libraries(${EXAMPLE_NAME} ${NAME})