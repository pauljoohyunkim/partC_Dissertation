cmake_minimum_required(VERSION 3.25.1)

project(gradient-flow)

add_library(SolverLib STATIC src/solver.cpp)

add_executable(L2 src/L2/L2.cpp)
target_link_libraries(L2 -lm SolverLib MathObjectsLib GeometricObjectsLib)

