cmake_minimum_required(VERSION 3.25.1)

project(CurveRepulsion LANGUAGES CUDA CXX)

add_library(GeometricObjects STATIC src/geometric-objects.cu)
add_library(Solver STATIC src/solver.cu src/geometric-objects.cu)

set(FOURIER_BIN "curve-repulsion-fourier")
set(CIRCLE_BIN "curve-repulsion-circle")

add_executable(${FOURIER_BIN} src/Curve-Repulsion/${FOURIER_BIN}.cu src/solver.cu src/geometric-objects.cu src/export.cu)
set_target_properties(${FOURIER_BIN} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
add_executable(${CIRCLE_BIN} src/Curve-Repulsion/${CIRCLE_BIN}.cu src/solver.cu src/geometric-objects.cu src/export.cu)
set_target_properties(${CIRCLE_BIN} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

#add_custom_target(upload1 DEPENDS curve-repulsion-1
#    COMMAND scp curve-repulsion-1 kim@gate.maths.ox.ac.uk:/home/kim/Documents/TEMP/curverepulsion1/)
#add_custom_target(upload2 DEPENDS curve-repulsion-1
#    COMMAND scp curve-repulsion-1 kim@gate.maths.ox.ac.uk:/home/kim/Documents/TEMP/curverepulsion2/)
#add_custom_target(upload3 DEPENDS curve-repulsion-1
#    COMMAND scp curve-repulsion-1 kim@gate.maths.ox.ac.uk:/home/kim/Documents/TEMP/curverepulsion3/)
