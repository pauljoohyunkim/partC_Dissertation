cmake_minimum_required(VERSION 3.21.1)

project(CurveRepulsionFunctionMinimize LANGUAGES CUDA CXX)

set(FOURIER_BIN "fourier-series-curve-repulsion-example-1")


add_executable(${FOURIER_BIN} src/FourierSeriesCurveRepulsion/${FOURIER_BIN}.cu src/curve.cu src/curverepulsion.cu)
set_target_properties(${FOURIER_BIN} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
