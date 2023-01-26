cmake_minimum_required(VERSION 3.25.1)

project(CurveRepulsion LANGUAGES CUDA CXX)

add_library(GeometricObjects STATIC ../src/geometric-objects.cu)
add_library(Solver STATIC ../src/solver.cu)

add_executable(cuCurveTest cuCurveTest.cu)
set_target_properties(cuCurveTest PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(cuCurveTest GeometricObjects)

add_executable(crossTest crossTest.cu ../src/solver.cu)
set_target_properties(crossTest PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
#target_link_libraries(crossTest Solver)
