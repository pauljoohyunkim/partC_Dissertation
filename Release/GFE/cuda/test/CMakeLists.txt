cmake_minimum_required(VERSION 3.25.1)

project(CurveRepulsion LANGUAGES CUDA CXX)


add_executable(loadTest loadTest.cu ../src/solver.cu)
set_target_properties(loadTest PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
