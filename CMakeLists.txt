cmake_minimum_required (VERSION 2.8)
project (gyper)

# The version number.
set (gyper_VERSION_MAJOR 0)
set (gyper_VERSION_MINOR 1)

# SeqAn options
set(SEQAN_BUILD_SYSTEM "SEQAN_RELEASE_LIBRARY" CACHE TYPE STRING)
set(CMAKE_BUILD_TYPE "Release" CACHE TYPE STRING)

# C++11
if(${GYPER_DEBUG})
  message("Debugging is enabled.")
  set (CMAKE_CXX_FLAGS "-g -W -Wall -Wfatal-errors -pedantic -Wno-variadic-macros -lrt -O2 -std=c++11 -DSEQAN_USE_HTSLIB=1 -DSEQAN_HAS_ZLIB=1 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=1")
else()
  message("Debugging is disabled.")
  set (CMAKE_CXX_FLAGS "-Wfatal-errors -lrt -std=c++11 -DSEQAN_USE_HTSLIB=1 -DSEQAN_HAS_ZLIB=1 -O3 -DNDEBUG -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0")
endif()

# Set module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/utilities/cmake_modules/")

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/src/constants.hpp.in"
  "${PROJECT_BINARY_DIR}/constants.hpp"
  )

include_directories("${PROJECT_BINARY_DIR}")

# Find external libraries
message (STATUS "Checking for zlib")
find_package(ZLIB 1.2.0 REQUIRED)
include_directories( ${ZLIB_INCLUDE_DIRS} )

message (STATUS "Checking for Boost")
find_package(Boost 1.57.0 REQUIRED)
include_directories( ${Boost_INCLUDE_DIRS} )

# Find SeqAn
# find_package(SeqAn REQUIRED)
# 
# if (NOT SEQAN_FOUND)
# 	add_subdirectory(utilities/seqan)
# 	find_package(SeqAn REQUIRED)
# endif(NOT SEQAN_FOUND)

# Make htslib
include(ExternalProject)
ExternalProject_Add(
  htslib
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/htslib
  BUILD_IN_SOURCE 1
  PREFIX ${CMAKE_BINARY_DIR}/htslib
  CONFIGURE_COMMAND ""
  BUILD_COMMAND $(MAKE) -C ${CMAKE_CURRENT_SOURCE_DIR}/htslib
  INSTALL_COMMAND ""
)

link_directories(${CMAKE_CURRENT_SOURCE_DIR}/htslib)
ExternalProject_Get_Property(htslib source_dir)
set(htslib_INCLUDE_DIRS ${source_dir} CACHE INTERNAL "Path to include folder for htslib.")
set(HTSlib_LIBRARY ${source_dir}/libhts.so CACHE INTERNAL "Path to htslib library files.")

# CMake SeqAn
add_subdirectory(seqan)
include_directories(seqan/include)

message("SeqAn include path is: ${SEQAN_INCLUDE_PATH}")
include_directories(${SEQAN_INCLUDE_DIRS})
include_directories(${SEQAN_INCLUDE_PATH})
# set(LIBS ${LIBS} ${SEQAN_LIBRARIES})

set(EXECUTABLE_OUTPUT_PATH ../bin)
add_subdirectory(src)

# Add vcflib
# add_subdirectory(vcflib)

message("Including HTSlib: " ${htslib_INCLUDE_DIRS})
include_directories(${htslib_INCLUDE_DIRS})

## Testing
include_directories(${CURRENT_SOURCE_DIR}/test)
# Include catch
include_directories(catch/include)
enable_testing(true)
add_subdirectory(test)

# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/docs/Doxyfile @ONLY)
add_custom_target(doc
  ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/docs/Doxyfile
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs
  COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif(DOXYGEN_FOUND)

