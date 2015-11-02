include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_executable(runUniformGridTests
        uniformGrid.cpp)

target_link_libraries(runUniformGridTests gtest gtest_main)
target_link_libraries(runUniformGridTests nufd)