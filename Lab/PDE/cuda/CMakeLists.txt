cmake_minimum_required(VERSION 3.25.1)

project(CurveRepulsion LANGUAGES CUDA CXX)

add_library(GeometricObjects STATIC src/geometric-objects.cu)
add_library(Solver STATIC src/solver.cu src/geometric-objects.cu)

add_executable(curve-repulsion-1 src/Curve-Repulsion/curve-repulsion-1.cu src/solver.cu src/geometric-objects.cu src/export.cu)
set_target_properties(curve-repulsion-1 PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(curve-repulsion-1 matplot)

