cmake_minimum_required(VERSION 2.6)
project(hf_solver)

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

add_library(hf_solver SHARED
	gto.cpp gto_eval.cpp orbitals.cpp system.cpp default.hpp hf.cpp basis_set.cpp
	common.hpp gto.hpp gto_eval.hpp math_helper.hpp orbitals.hpp system.hpp hf.hpp basis_set.hpp
)

#target_link_libraries(hf_solver
#	gsl gslcblas
#)
enable_testing()
add_subdirectory(example)
add_subdirectory(test)
