cmake_minimum_required(VERSION 3.12)
project(simulations)

set(CMAKE_CXX_STANDARD 11)

add_executable(simulations main.cpp src/SimulationEngine.cpp src/SimulationEngine.h)