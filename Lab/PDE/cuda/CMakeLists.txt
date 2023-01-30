cmake_minimum_required(VERSION 3.25.1)

project(CurveRepulsion LANGUAGES CUDA CXX)

add_library(GeometricObjects STATIC src/geometric-objects.cu)
add_library(Solver STATIC src/solver.cu src/geometric-objects.cu)

add_executable(curve-repulsion-1 src/Curve-Repulsion/curve-repulsion-1.cu src/solver.cu src/geometric-objects.cu src/export.cu)
set_target_properties(curve-repulsion-1 PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

add_custom_target(upload1 DEPENDS curve-repulsion-1
    COMMAND scp curve-repulsion-1 kim@gate.maths.ox.ac.uk:/home/kim/Documents/TEMP/curverepulsion1/)
add_custom_target(upload2 DEPENDS curve-repulsion-1
    COMMAND scp curve-repulsion-1 kim@gate.maths.ox.ac.uk:/home/kim/Documents/TEMP/curverepulsion2/)
add_custom_target(upload3 DEPENDS curve-repulsion-1
    COMMAND scp curve-repulsion-1 kim@gate.maths.ox.ac.uk:/home/kim/Documents/TEMP/curverepulsion3/)
