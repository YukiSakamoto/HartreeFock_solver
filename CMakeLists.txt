cmake_minimum_required(VERSION 2.6)
project(hf_solver)

add_library(hf_solver SHARED
	gto.cpp gto_eval.cpp orbitals.cpp system.cpp hf.cpp
	common.hpp gto.hpp gto_eval.hpp math_helper.hpp orbitals.hpp system.hpp hf.hpp
)

#target_link_libraries(hf_solver
#	gsl gslcblas
#)
enable_testing()
add_subdirectory(example)
add_subdirectory(test)
