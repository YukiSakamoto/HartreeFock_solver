cmake_minimum_required(VERSION 2.6)
project(run_rhf)

include_directories(..)
add_executable(run_rhf main.cpp)

target_link_libraries(run_rhf
	hf_solver
)
