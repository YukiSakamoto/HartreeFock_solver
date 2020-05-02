cmake_minimum_required(VERSION 2.6)
project(hf_solver)

add_executable(MOSolver 
	main.cpp)

target_link_libraries(MOSolver
	hf_solver
)

#add_subdirectory(example)
#add_subdirectory(test)
