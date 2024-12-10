#
# @file      CMakeLists.txt
#
# @author    David Bayer \n
#            Faculty of Information Technology \n
#            Brno University of Technology \n
#            ibayer@fit.vutbr.cz
#
# @brief     PCG Assignment 2
#
# @version   2023
#
# @date      04 October   2023, 09:00 (created) \n
#

cmake_minimum_required(VERSION 3.8)

if(NOT DEFINED ACC_TARGET)
  set(ACC_TARGET "gpu")
endif()

message(STATUS "ACC_TARGET: ${ACC_TARGET}")

project(Assignment1 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD               17)
set(CMAKE_CXX_STANDARD_REQUIRED      ON)

if(NOT CMAKE_CXX_COMPILER_ID MATCHES "NVHPC")
  message(FATAL_ERROR "Unsupported compiler. Use NVHPC compiler.")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -acc=${ACC_TARGET} -Minfo=accel -mp=multicore")
if(ACC_TARGET MATCHES "gpu")
  set(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG} -gpu=ccall,debug,pinned")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -gpu=ccall,fastmath,pinned")
endif()

find_package(NVHPC REQUIRED COMPONENTS HOSTUTILS)
link_directories(NVHPC::HOSTUTILS)

find_package(HDF5 REQUIRED COMPONENTS C)

include_directories(Commons)
link_libraries(HDF5::HDF5)

# Gen
add_executable(gen Commons/gen.cpp)

# Cpu
add_executable(nbodyCpu Cpu/main.cpp Cpu/nbody.cpp Commons/h5Helper.cpp)

# Step0
add_executable(nbody0 Step0/main.cpp Step0/nbody.cpp Commons/h5Helper.cpp)

# Step1
add_executable(nbody1 Step1/main.cpp Step1/nbody.cpp Commons/h5Helper.cpp)

# Step2
add_executable(nbody2 Step2/main.cpp Step2/nbody.cpp Commons/h5Helper.cpp)

# Step3
add_executable(nbody3 Step3/main.cpp Step3/nbody.cpp Commons/h5Helper.cpp)
