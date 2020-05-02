


add_executable(MOSolver 
	main.cpp)

target_link_libraries(MOSolver
	hf_solver
)

file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/examples)
file(GLOB EXAMPLES ${CMAKE_CURRENT_SOURCE_DIR}/example/*)
file(COPY ${EXAMPLES} DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/examples)
