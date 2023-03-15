cmake_minimum_required(VERSION 3.25.1)

project(CurveRepulsion LANGUAGES CUDA CXX)

set(FOURIER_BIN "curve-repulsion-fourier")

add_executable(${FOURIER_BIN} src/Curve-Repulsion/${FOURIER_BIN}.cu src/solver.cu src/export.cu)
set_target_properties(${FOURIER_BIN} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
