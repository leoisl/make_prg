cmake_minimum_required(VERSION 3.3)
project(make_prg)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")



################################################################################
# getting git version
################################################################################
#from http://stackoverflow.com/questions/1435953/how-can-i-pass-git-sha1-to-compiler-as-definition-using-cmake
exec_program(
        "git"
        ${CMAKE_CURRENT_SOURCE_DIR}
        ARGS "rev-parse --short HEAD"
        OUTPUT_VARIABLE VERSION_SHA1 )

add_definitions( -DGIT_SHA1="${VERSION_SHA1}" )
################################################################################
# getting git version
################################################################################



################################################################################
# Boost configuration
################################################################################
set(Boost_USE_STATIC_LIBS ON)
find_package(Boost 1.62.0 COMPONENTS system filesystem iostreams log thread date_time REQUIRED)
link_directories(${Boost_LIBRARY_DIR})
################################################################################
# Boost configuration
################################################################################



################################################################################
# GATB configuration
################################################################################
# The following lines are inspired by https://github.com/GATB/bcalm/blob/master/CMakeLists.txt
SET(GATB_CORE_HOME  ${PROJECT_SOURCE_DIR}/gatb-core)

################################################################################
# Define cmake modules directory
################################################################################
FOREACH (path "${CMAKE_CURRENT_SOURCE_DIR}/cmake" "${GATB_CORE_HOME}/gatb-core/cmake")
    IF (EXISTS "${path}")
        SET (CMAKE_MODULE_PATH  "${CMAKE_MODULE_PATH}" "${path}")
    ENDIF()
ENDFOREACH(path)


################################################################################
# SUPPORTED KMER SIZES
################################################################################
# One can uncomment this line and set the wanted values
# TODO: if you change here, search the project for KSIZE_LIST and change the
# code accordingly
set (KSIZE_LIST "32   64   96  128")


################################################################################
# THIRD PARTIES
################################################################################

# We don't want to install some GATB-CORE artifacts
SET (GATB_CORE_EXCLUDE_TOOLS     1) # no need to compile dbgh5, etc..
SET (GATB_CORE_EXCLUDE_TESTS     1)
SET (GATB_CORE_EXCLUDE_EXAMPLES  1)

# GATB CORE
include (GatbCore)

# We add the compilation options for the library
#add_definitions (${gatb-core-flags}) #we don't want to add GATB's compilations options back here

# We add the gatb-core include directory
#include_directories (${gatb-core-includes}) #if we add this, then the boost in gatb-core/gatb-core/thirdparty conflicts with pre-installed boost and compilation fails
#just add what we need without GATB boost
SET (GATB_SOURCE_DIR "${PROJECT_SOURCE_DIR}/gatb-core/gatb-core")
SET (GATB_BINARY_DIR "${PROJECT_BINARY_DIR}/ext/gatb-core")
SET (gatb-core-includes ${GATB_BINARY_DIR}/include  ${GATB_BINARY_DIR}/include/${CMAKE_BUILD_TYPE} ${GATB_SOURCE_DIR}/src ${gatb-core-extra-libraries-inc})
include_directories(${gatb-core-includes})

# We add the path for extra libraries
link_directories (${gatb-core-extra-libraries-path})
################################################################################
# GATB configuration
################################################################################




################################################################################
# Building the target
################################################################################
include_directories(${PROJECT_SOURCE_DIR} "${PROJECT_SOURCE_DIR}/BooM/include")
file (GLOB_RECURSE  SOURCE_FILES  ${PROJECT_SOURCE_DIR}/cpp/*.cpp)
add_executable(make_prg ${SOURCE_FILES})
target_link_libraries(make_prg ${Boost_LIBRARIES} ${gatb-core-libraries} -static-libgcc -static-libstdc++)
################################################################################
# Building the target
################################################################################