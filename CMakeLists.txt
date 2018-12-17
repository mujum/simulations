cmake_minimum_required(VERSION 3.12)
project(simulations)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_COMPILER g++)

find_package(HDF5 REQUIRED COMPONENTS C CXX)

include_directories(${HDF5_INCLUDE_DIRS})

add_executable(simulations main.cpp
        src/SimulationEngine.cpp src/SimulationEngine.h
        src/Visualization.cpp src/Visualization.h
        )

target_link_libraries( simulations
        ${HDF5_LIBRARIES}
        )