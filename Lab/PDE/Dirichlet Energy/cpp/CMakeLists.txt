cmake_minimum_required(VERSION 3.25.1)

project(gradient-flow)

add_executable(L2 src/L2/L2.cpp src/L2/solvers.cpp)
target_link_libraries(L2 -lm -lmatplot)
