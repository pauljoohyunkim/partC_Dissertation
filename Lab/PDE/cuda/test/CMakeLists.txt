cmake_minimum_required(VERSION 3.25.1)

project(CurveRepulsion LANGUAGES CUDA CXX)

add_library(GeometricObjects STATIC ../src/geometric-objects.cu)

add_executable(cuCurveTest cuCurveTest.cu)
set_target_properties(cuCurveTest PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(cuCurveTest GeometricObjects)

