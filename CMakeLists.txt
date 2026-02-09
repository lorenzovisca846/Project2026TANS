cmake_minimum_required(VERSION 3.10)
project(Project2026TANS)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})

include_directories(${CMAKE_SOURCE_DIR}/classes)

add_subdirectory(classes)

add_executable(simulation main_sim.cpp)
target_link_libraries(simulation ProjectClasses ${ROOT_LIBRARIES})