cmake_minimum_required (VERSION 2.8)
project (gyper)

# The version number.
set (gyper_VERSION_MAJOR 0)
set (gyper_VERSION_MINOR 1)

# C++14
set (CMAKE_CXX_FLAGS "-O3 -std=c++14 -DSEQAN_HAS_ZLIB=1 -DNDEBUG -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0")

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
include_directories(src)

# Find files
file(GLOB_RECURSE gyper_SOURCES "src/*.cpp")
add_executable(gyper ${gyper_SOURCES})

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

# Installs gyper to /usr/local/bin using [sudo] make install
install (TARGETS gyper DESTINATION bin)
