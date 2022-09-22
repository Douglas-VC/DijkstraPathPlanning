cmake_minimum_required(VERSION 3.5.1)
project(DijkstraPathPlanning)
set(CMAKE_CXX_STANDARD 11)

find_package(PCL 1.7 REQUIRED COMPONENTS common io)

include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

SET(TARGET_H
        coder_array.h
        Dijkstra.h
        rtwtypes.h
        find.h
        rt_nonfinite.h
        tmwtypes.h
        )

# Local source files here
SET(TARGET_SRC
        Dijkstra.cpp
        find.cpp
        rt_nonfinite.cpp
        )

add_executable(DijkstraPathPlanning DijkstraPathPlanning.cpp ${TARGET_H} ${TARGET_SRC})

target_link_libraries(DijkstraPathPlanning ${PCL_LIBRARIES})
