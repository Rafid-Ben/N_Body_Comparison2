cmake_minimum_required(VERSION 3.0.0)
project(NBodySim VERSION 0.1.0)
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -std=c++11)

include(CTest)
enable_testing()
#set(CMAKE_SYSTEM_PROCESSOR arm64)
FILE(GLOB HEADER_FILES *.h)
FILE(GLOB SOURCES *c *.cpp)
add_executable(NBodySim_test test_cases.cpp ${HEADER_FILES} ${SOURCES})
set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)


#set (CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_STANDARD 25)
#set(CMAKE_CXX_STANDARD_REQUIRED on)

