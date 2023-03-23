cmake_minimum_required(VERSION 3.25.1)

project(CurveRepulsion LANGUAGES CUDA CXX)

#set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS}; -gencode arch=compute_86,code=sm_86")
#set(CMAKE_VERBOSE_MAKEFILE TRUE)
set(CUDA_GENCODE_NUM 86)


add_executable(loadTest loadTest.cu ../src/solver.cu)
set_target_properties(loadTest PROPERTIES CUDA_SEPARABLE_COMPILATION ON CUDA_ARCHITECTURES ${CUDA_GENCODE_NUM})

add_executable(addTest addTest.cu ../src/solver.cu)
set_target_properties(addTest PROPERTIES CUDA_SEPARABLE_COMPILATION ON CUDA_ARCHITECTURES ${CUDA_GENCODE_NUM})

add_executable(accessTest accessTest.cu ../src/solver.cu)
set_target_properties(accessTest PROPERTIES CUDA_SEPARABLE_COMPILATION ON CUDA_ARCHITECTURES ${CUDA_GENCODE_NUM})

add_executable(vectorTest vectorTest.cu ../src/solver.cu ../src/vector.cu)
set_target_properties(vectorTest PROPERTIES CUDA_SEPARABLE_COMPILATION ON CUDA_ARCHITECTURES ${CUDA_GENCODE_NUM})

add_executable(vectorFromTensorTest vectorFromTensorTest.cu ../src/solver.cu ../src/vector.cu)
set_target_properties(vectorFromTensorTest PROPERTIES CUDA_SEPARABLE_COMPILATION ON CUDA_ARCHITECTURES ${CUDA_GENCODE_NUM})
