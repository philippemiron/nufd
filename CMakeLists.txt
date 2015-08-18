cmake_minimum_required(VERSION 3.2)
project(nufd)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    nufd.cpp
    nufd.h
    )

add_subdirectory(tests/gtest-1.7.0/)
include_directories(tests/gtest-1.7.0/include/)

# example executable
add_executable(nufd ${SOURCE_FILES} example.cpp)

# gtest executable
add_executable(unitTest ${SOURCE_FILES} tests/uniformGrid.cpp)
target_link_libraries(unitTest gtest gtest_main)