cmake_minimum_required(VERSION 3.25.1)

project(CurveRepulsion LANGUAGES CUDA CXX)

set_target_properties(curve-repulsion-1 PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
#add_library(SolverLib STATIC src/solver.cpp)

add_executable(curve-repulsion-1 src/Curve-Repulsion/curve-repulsion-1.cu)
#target_link_libraries(L2 -lm SolverLib MathObjectsLib GeometricObjectsLib)

