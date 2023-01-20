cmake_minimum_required(VERSION 3.25.1)

project(gradient-flow)

add_library(SolverLib STATIC src/solver.cpp)

add_executable(L2 src/L2/L2.cpp)
add_executable(H1 src/H1/H1.cpp)
add_executable(H-1 src/H-1/H-1.cpp)
add_executable(porous src/Porous-Medium/porous.cpp)
target_link_libraries(L2 -lm SolverLib)
target_link_libraries(H1 -lm SolverLib)
target_link_libraries(H-1 -lm SolverLib)
target_link_libraries(porous -lm SolverLib)

