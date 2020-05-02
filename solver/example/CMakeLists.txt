cmake_minimum_required(VERSION 2.6)
project(example_binary)

include_directories(..)
add_executable(example_HeH main.cpp)

target_link_libraries(example_HeH
	hf_solver
)

add_executable(example_CH4 ch4.cpp)
target_link_libraries(example_CH4 hf_solver)

file(GLOB RESOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/example_resources/*.dat)
file(COPY ${RESOURCE_FILES} DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

