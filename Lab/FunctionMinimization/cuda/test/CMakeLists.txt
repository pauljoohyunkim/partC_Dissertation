cmake_minimum_required(VERSION 3.21.1)

project(CurveRepulsionTest LANGUAGES CUDA CXX)

add_library(fouriercurve STATIC ../src/curve.cu)

add_executable(cuCurveTest cuCurveTest.cu ../src/curve.cu)
set_target_properties(cuCurveTest PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

