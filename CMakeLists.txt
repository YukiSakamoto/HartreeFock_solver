cmake_minimum_required(VERSION 2.6)
project(hf_solver CXX)

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -march=native")
set(CMAKE_CXX_FLAGS_MINSIZEREL "-Os -s -march=native")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

enable_testing()
add_subdirectory(solver)
add_subdirectory(driver)

