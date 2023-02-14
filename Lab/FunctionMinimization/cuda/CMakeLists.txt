cmake_minimum_required(VERSION 3.21.1)

project(CurveRepulsionFunctionMinimize LANGUAGES CUDA CXX)

set(FOURIER_BIN_1 "fourier-series-curve-repulsion-example-1")
set(FOURIER_BIN_2 "fourier-series-curve-repulsion-example-2")


add_executable(${FOURIER_BIN_1} src/FourierSeriesCurveRepulsion/${FOURIER_BIN_1}.cu src/curve.cu src/curverepulsion.cu src/export.cu)
set_target_properties(${FOURIER_BIN_1} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

add_executable(${FOURIER_BIN_2} src/FourierSeriesCurveRepulsion/${FOURIER_BIN_2}.cu src/curve.cu src/curverepulsion.cu src/export.cu)
set_target_properties(${FOURIER_BIN_2} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
