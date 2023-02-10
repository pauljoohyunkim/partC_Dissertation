cmake_minimum_required(VERSION 3.25.1)

project(CurveRepulsionFunctionMinimize LANGUAGES CUDA CXX)

add_library(GeometricObjects STATIC src/geometric-objects.cu)
add_library(Solver STATIC src/solver.cu src/geometric-objects.cu)

set(FOURIER_BIN "curve-repulsion-fourier")


add_executable(${FOURIER_BIN} src/Curve-Repulsion/${FOURIER_BIN}.cu src/solver.cu)
set_target_properties(${FOURIER_BIN} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
