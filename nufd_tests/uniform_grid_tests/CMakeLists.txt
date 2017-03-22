include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_executable(UniformGridTests
        uniformGrid.cpp)

target_link_libraries(UniformGridTests gtest gtest_main)
target_link_libraries(UniformGridTests nufd)