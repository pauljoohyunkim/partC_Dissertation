cmake_minimum_required(VERSION 3.25.1)

project(CurveRepulsion LANGUAGES CUDA CXX)

set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS}; -gencode arch=compute_61,code=sm_61")
set(CMAKE_VERBOSE_MAKEFILE TRUE)


add_executable(loadTest loadTest.cu ../src/solver.cu)
set_target_properties(loadTest PROPERTIES CUDA_SEPARABLE_COMPILATION ON CUDA_ARCHITECTURES 61)

add_executable(addTest addTest.cu ../src/solver.cu)
set_target_properties(addTest PROPERTIES CUDA_SEPARABLE_COMPILATION ON CUDA_ARCHITECTURES 61)
