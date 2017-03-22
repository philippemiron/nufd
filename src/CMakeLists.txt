add_definitions(-std=c++11)

set(HEADER_FILES nufd.h)

set(SOURCE_FILES nufd.cpp)

add_library(nufd STATIC ${SOURCE_FILES} ${HEADER_FILES})
