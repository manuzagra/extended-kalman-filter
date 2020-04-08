project(Extended_Kalman_Filter)

cmake_minimum_required (VERSION 3.5)

add_definitions(-std=c++11)

set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS "${CXX_FLAGS}")

set(sources
	src/main.cpp
	src/tools.cpp
	src/FusionEKF.cpp
	src/kalman_filter.cpp)

add_executable(ExtendedKF ${sources})

target_link_libraries(ExtendedKF z ssl uv uWS)
