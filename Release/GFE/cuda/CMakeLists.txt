cmake_minimum_required(VERSION 3.25.1)

project(CurveRepulsion LANGUAGES CUDA CXX)

set(CUDA_GENCODE_NUM 86)
set(FIGURE_8 "figure8")

add_executable(${FIGURE_8} src/demos/${FIGURE_8}.cu src/solver.cu src/export.cu src/vector.cu src/tpe.cu src/export.cu)
set_target_properties(${FIGURE_8} PROPERTIES CUDA_SEPARABLE_COMPILATION ON CUDA_ARCHITECTURES ${CUDA_GENCODE_NUM})
