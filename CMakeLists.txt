cmake_minimum_required(VERSION 3.2)
project(Nufd)

add_definitions(-std=c++11)


include_directories(src)
set(SOURCE_FILES example.cpp)

# example executable
add_executable(nufd_example ${SOURCE_FILES})

# add cmake subdirectories
add_subdirectory(src)
add_subdirectory(nufd_tests)

target_link_libraries(nufd_example nufd)