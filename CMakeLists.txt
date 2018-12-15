cmake_minimum_required(VERSION 3.12)
project(simulations)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_COMPILER g++)

#include_directories(${HDF5_INCLUDE_DIRS})

#find_package(HDF5 REQUIRED COMPONENTS C CXX)
find_package(HDF5 REQUIRED COMPONENTS C CXX)

#include_directories(/disk/learnspace/_laboratory/bin/hdf5-1.10.4-linux-centos7-x86_64-gcc485-shared/include)
#link_directories(/disk/learnspace/_laboratory/bin/hdf5-1.10.4-linux-centos7-x86_64-gcc485-shared/lib)

include_directories(${HDF5_INCLUDE_DIRS})

add_executable(simulations main.cpp src/SimulationEngine.cpp src/SimulationEngine.h)

target_link_libraries( simulations ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES} )# libhdf5.so libhdf5_cpp.so libhdf5_hl.so libhdf5_hl_cpp.so libz.so libsz.so)