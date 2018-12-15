cmake_minimum_required(VERSION 3.12)
project(simulations)

set(CMAKE_CXX_STANDARD 11)


include_directories(/disk/learnspace/_laboratory/bin/hdf5-1.10.4-linux-centos7-x86_64-gcc485-shared/include)

link_directories(/disk/learnspace/_laboratory/bin/hdf5-1.10.4-linux-centos7-x86_64-gcc485-shared/lib)

add_executable(simulations main.cpp src/SimulationEngine.cpp src/SimulationEngine.h)