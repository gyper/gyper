cmake_minimum_required (VERSION 2.8.8)
project (gyper)

# The version number.
set (gyper_VERSION_MAJOR 0)
set (gyper_VERSION_MINOR 1)

# C++14
set (CMAKE_CXX_FLAGS "-W -Wall -Wno-long-long -pedantic -Wno-variadic-macros -lrt -O3 -std=c++14 -DSEQAN_HAS_ZLIB=1 -DNDEBUG -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0")

# Set module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/src/constants.hpp.in"
  "${PROJECT_BINARY_DIR}/constants.hpp"
  )

include_directories("${PROJECT_BINARY_DIR}")

# Binary directory
set(EXT_PROJECTS_DIR ${PROJECT_SOURCE_DIR}/test)

# include_directories(src)
# add_subdirectory(test/catch)
# include_directories(${CATCH_INCLUDE_DIR} src test)
# enable_testing(true)  # Enables unit-testing.

# Define source and test files
set(GYPER_SOURCE_FILES src/graph.cpp src/graph_align.cpp src/graph_builder.cpp src/graph_io.cpp src/graph_kmerify.cpp)
set(GYPER_TEST_FILES test_gyper test/test.cpp test/test_graph_builder.cpp test/test_graph_align.cpp)

# Add executables
add_executable(gyper src/main.cpp src/graph.cpp src/graph_align.cpp src/graph_builder.cpp src/graph_io.cpp src/graph_kmerify.cpp)
add_executable(test_gyper test/test.cpp test/test_graph_builder.cpp test/test_graph_align.cpp src/graph.cpp src/graph_align.cpp src/graph_builder.cpp src/graph_io.cpp src/graph_kmerify.cpp)

# Find libraries
message (STATUS "Checking for zlib")
find_package(ZLIB REQUIRED)
include_directories( ${ZLIB_INCLUDE_DIRS} )
target_link_libraries( gyper ${ZLIB_LIBRARIES} )

message (STATUS "Checking for Boost")
find_package(Boost 1.56.0 REQUIRED)
include_directories( ${Boost_INCLUDE_DIRS} )
target_link_libraries( gyper ${Boost_LIBRARIES} )

find_package(SeqAn REQUIRED)
include_directories(${SEQAN_INCLUDE_DIRS})
message (STATUS "Added include directories: ${SEQAN_INCLUDE_DIRS}")
set(LIBS ${LIBS} ${SEQAN_LIBRARIES})

# Add SeqAn definitions set by find_package.
add_definitions (${SEQAN_DEFINITIONS})

if (NOT SEQAN_FOUND)
	message (ERROR "SeqAn not found!")
endif(NOT SEQAN_FOUND)

# Add the executable
set(EXECUTABLE_OUTPUT_PATH bin)
target_link_libraries(gyper ${LIBS})
target_link_libraries(test_gyper ${LIBS})

# Testing
include_directories(${EXT_PROJECTS_DIR})
add_subdirectory(${EXT_PROJECTS_DIR}/catch)
include_directories(${CATCH_INCLUDE_DIR} ${COMMON_INCLUDES})
enable_testing(true)

# Installs gyper to /usr/local/bin using [sudo] make install
install (TARGETS gyper DESTINATION bin)
