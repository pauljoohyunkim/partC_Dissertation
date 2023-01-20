cmake_minimum_required(VERSION 3.25.1)

project(tests)

add_library(MathObjectsLib STATIC ../src/math-objects.cpp)
add_executable(matrix-test matrix-test.cpp)
target_link_libraries(matrix-test MathObjectsLib)
