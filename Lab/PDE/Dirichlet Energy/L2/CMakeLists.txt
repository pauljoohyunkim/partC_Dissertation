cmake_minimum_required(VERSION 3.25.1)

project(L2)

add_executable(L2 src/L2.cpp src/solvers.cpp)
target_link_libraries(L2 -lm -lmatplot)
