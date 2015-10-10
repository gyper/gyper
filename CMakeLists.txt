cmake_minimum_required (VERSION 2.8)
project (gyper)

# The version number.
set (gyper_VERSION_MAJOR 0)
set (gyper_VERSION_MINOR 1)

# C++11
set (CMAKE_CXX_FLAGS "-W -Wall -Wno-long-long -pedantic -Wno-variadic-macros -lrt -O3 -std=c++11 -DSEQAN_HAS_ZLIB=1 -DNDEBUG -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0")

# Set module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/src/constants.hpp.in"
  "${PROJECT_BINARY_DIR}/constants.hpp"
  )

include_directories("${PROJECT_BINARY_DIR}")

# Find external libraries
message (STATUS "Checking for zlib")
find_package(ZLIB REQUIRED)
include_directories( ${ZLIB_INCLUDE_DIRS} )

message (STATUS "Checking for Boost")
find_package(Boost 1.57.0 REQUIRED)
include_directories( ${Boost_INCLUDE_DIRS} )

# Find SeqAn
find_package(SeqAn REQUIRED)

if (NOT SEQAN_FOUND)
	add_subdirectory(cmake/seqan)
	find_package(SeqAn REQUIRED)
endif(NOT SEQAN_FOUND)

message("SEQAN_INCLUDE_PATH is ${SEQAN_INCLUDE_PATH}")
include_directories(${SEQAN_INCLUDE_DIRS})
include_directories(${SEQAN_INCLUDE_PATH})
# set(LIBS ${LIBS} ${SEQAN_LIBRARIES})

set(EXECUTABLE_OUTPUT_PATH ../bin)
add_subdirectory(src)

# Testing
set(EXT_PROJECTS_DIR ${PROJECT_SOURCE_DIR}/test)
include_directories(${EXT_PROJECTS_DIR})
add_subdirectory(${EXT_PROJECTS_DIR}/catch)
include_directories(${CATCH_INCLUDE_DIR} ${COMMON_INCLUDES})
enable_testing(true)

add_subdirectory(test)
