cmake_minimum_required(VERSION 3.25.1)

project(gradient-flow)

add_library(SolverLib STATIC src/solver.cpp)
add_library(MathObjectsLib STATIC src/math-objects.cpp)
add_library(GeometricObjectsLib STATIC src/geometric-objects.cpp)
add_library(DiscreteKernelLib STATIC src/Curve-Repulsion/discrete-kernel.cpp)

add_executable(L2 src/L2/L2.cpp)
add_executable(H1 src/H1/H1.cpp)
add_executable(H-1 src/H-1/H-1.cpp)
add_executable(porous src/Porous-Medium/porous.cpp)
add_executable(curve-repulsion-1 src/Curve-Repulsion/curve-repulsion-1.cpp)
target_link_libraries(L2 -lm SolverLib)
target_link_libraries(H1 -lm SolverLib)
target_link_libraries(H-1 -lm SolverLib)
target_link_libraries(porous -lm SolverLib)
target_link_libraries(curve-repulsion-1 -lm DiscreteKernelLib MathObjectsLib GeometricObjectsLib -lmatplot)

