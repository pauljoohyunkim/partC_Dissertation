cmake_minimum_required(VERSION 3.25.1)

project(gradient-flow)

add_library(SolverLib STATIC src/solver.cpp)
add_library(MathObjectsLib STATIC src/math-objects.cpp)
add_library(GeometricObjectsLib STATIC src/geometric-objects.cpp)
add_library(ExportLib STATIC src/export.cpp)

set(CIRCLE_BIN "curve-repulsion-circle")
set(DIAMOND_BIN "curve-repulsion-diamond")

add_executable(L2 src/L2/L2.cpp)
add_executable(H1 src/H1/H1.cpp)
add_executable(H-1 src/H-1/H-1.cpp)
add_executable(porous src/Porous-Medium/porous.cpp)
add_executable(${CIRCLE_BIN} src/Curve-Repulsion/${CIRCLE_BIN}.cpp)
add_executable(${DIAMOND_BIN} src/Curve-Repulsion/${DIAMOND_BIN}.cpp)
target_link_libraries(L2 -lm SolverLib MathObjectsLib GeometricObjectsLib)
target_link_libraries(H1 -lm SolverLib MathObjectsLib GeometricObjectsLib)
target_link_libraries(H-1 -lm SolverLib MathObjectsLib GeometricObjectsLib)
target_link_libraries(porous -lm SolverLib MathObjectsLib GeometricObjectsLib)
target_link_libraries(${CIRCLE_BIN} -lm SolverLib MathObjectsLib GeometricObjectsLib ExportLib)
target_link_libraries(${DIAMOND_BIN} -lm SolverLib MathObjectsLib GeometricObjectsLib ExportLib)

