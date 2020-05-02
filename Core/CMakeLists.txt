cmake_minimum_required(VERSION 2.6)
#project(hf_solver)

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

add_library(MOSolverCore SHARED
	GaussBasis.hpp System.hpp Core.hpp 
	GaussBasis.cpp
)

#enable_testing()
add_subdirectory(example)
#add_subdirectory(test)
