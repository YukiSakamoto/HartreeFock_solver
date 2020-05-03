cmake_minimum_required(VERSION 2.6)

add_library(hf_solver SHARED
	gto.cpp gto_eval.cpp orbitals.cpp system.cpp default.hpp hf.cpp basis_set.cpp
	common.hpp gto.hpp gto_eval.hpp math_helper.hpp orbitals.hpp system.hpp hf.hpp basis_set.hpp
)

option(USE_LOOP_OPT "Loop optimization code will be used" OFF)
if (USE_LOOP_OPT)
	add_definitions(-DLOOP_OPT)
endif()

# Copy the example files to output directory
set(CMAKE_CURRENT_BINARY_DIR ${CMAKE_BINARY_DIR})

add_subdirectory(example)
add_subdirectory(test)

