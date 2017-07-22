cmake_minimum_required(VERSION 2.6)
project(run_rhf)

include_directories(..)
add_executable(run_rhf main.cpp)

target_link_libraries(run_rhf
	hf_solver
)

add_executable(ch4 ch4.cpp)
target_link_libraries(ch4 hf_solver)

add_executable(ch4_dft ch4_dft.cpp)
target_link_libraries(ch4_dft hf_solver)
