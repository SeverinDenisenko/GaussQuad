cmake_minimum_required(VERSION 3.22)
project(numerical_task_8)

set(CMAKE_CXX_STANDARD 20)
add_compile_options("-Wall -Wextra -Wfloat-equal -Wundef -Wcast-align -Wwrite-strings -Wlogical-op -Wmissing-declarations -Wredundant-decls -Wshadow -Woverloaded-virtual -pedantic")

add_executable(solver gauss_solver/main.cpp gauss_solver/gauss.cpp gauss_solver/gauss.h)
add_executable(integrator gauss_integrator/main.cpp)
