cmake_minimum_required(VERSION 3.3)
project(make_prg)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    cpp/BuildPRG.cpp
    cpp/BuildPRG.h
    cpp/main.cpp
    cpp/PRGBoostGraph.h
    cpp/PRGBoostGraph.cpp)

add_executable(make_prg ${SOURCE_FILES})