cmake_minimum_required(VERSION 3.10)
project(SimulationProject)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(ROOT REQUIRED COMPONENTS Core Hist RIO Net Gpad Tree MathCore Physics)
include(${ROOT_USE_FILE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/simclasses)

set(SIM_SOURCES
    simclasses/Cylinder.cxx
    simclasses/MyPoint.cxx
    simclasses/Particle.cxx
    simclasses/SimRandom.cxx
)

root_generate_dictionary(G__Simulation 
    simclasses/Cylinder.h 
    simclasses/MyPoint.h 
    simclasses/Particle.h 
    simclasses/SimRandom.h 
    LINKDEF LinkDef.h
    OPTIONS -I${CMAKE_CURRENT_SOURCE_DIR}
)

add_library(SimClasses SHARED ${SIM_SOURCES} G__Simulation.cxx)
target_link_libraries(SimClasses PUBLIC ${ROOT_LIBRARIES})

add_executable(simulation main.cpp)

target_link_libraries(simulation SimClasses ${ROOT_LIBRARIES})