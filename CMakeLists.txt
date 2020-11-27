cmake_minimum_required(VERSION 3.17)
project(program_optimization)

set(CMAKE_CXX_STANDARD 20)

add_executable(program_optimization loop_unrolling.cpp polynomial_evaluation.cpp floating_point_error_mitigation.cpp
        my_solver.h main.cpp)