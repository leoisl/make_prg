cmake_minimum_required(VERSION 3.3)
project(make_prg)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

include_directories("${PROJECT_SOURCE_DIR}")

set(Boost_USE_STATIC_LIBS ON)
find_package(Boost 1.62.0 COMPONENTS system filesystem iostreams log thread date_time REQUIRED)
link_directories(${Boost_LIBRARY_DIR})

file (GLOB_RECURSE  SOURCE_FILES  ${PROJECT_SOURCE_DIR}/cpp/*.cpp)

add_executable(make_prg ${SOURCE_FILES})
target_link_libraries(make_prg ${Boost_LIBRARIES} # -static-libgcc -static-libstdc++
        )