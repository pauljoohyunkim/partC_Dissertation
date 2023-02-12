cmake_minimum_required(VERSION 3.21.1)

project(CurveRepulsionTest LANGUAGES CUDA CXX)

add_library(fouriercurve STATIC ../src/curve.cu)

add_executable(cuCurveTest cuCurveTest.cu ../src/curve.cu ../src/curverepulsion.cu)
set_target_properties(cuCurveTest PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

add_executable(cuCrossTest cuCrossTest.cu ../src/curve.cu)
set_target_properties(cuCrossTest PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

add_executable(exportTest exportTest.cu ../src/export.cu)
set_target_properties(exportTest PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
