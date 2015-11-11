cmake_minimum_required (VERSION 2.8)
project (gyper)

include(ExternalProject)

# The version number.
set (gyper_VERSION_MAJOR 0)
set (gyper_VERSION_MINOR 1)

# SeqAn options
set(SEQAN_BUILD_SYSTEM "SEQAN_RELEASE_LIBRARY" CACHE TYPE STRING)
set(CMAKE_BUILD_TYPE "Release" CACHE TYPE STRING)

# C++11
if(${GYPER_DEBUG})
  message(STATUS "Debugging is enabled.")
  set (CMAKE_CXX_FLAGS "-g -W -Wall -Wfatal-errors -pedantic -Wno-variadic-macros -lrt -O2 -std=c++11 -DSEQAN_USE_HTSLIB=1 -DSEQAN_HAS_ZLIB=1 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=1")
else()
  message(STATUS "Debugging is disabled.")
  set (CMAKE_CXX_FLAGS "-Wfatal-errors -lrt -std=c++11 -DSEQAN_USE_HTSLIB=1 -DSEQAN_HAS_ZLIB=1 -O3 -DNDEBUG -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0")
endif()

# Set module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/utilities/cmake_modules/")

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/include/gyper/constants.hpp.in"
  "${PROJECT_BINARY_DIR}/gyper/constants.hpp"
  )

include_directories("${PROJECT_BINARY_DIR}")

# Find external libraries
# zlib
message (STATUS "Checking for zlib")
find_package(ZLIB 1.2.0 REQUIRED)
include_directories( ${ZLIB_INCLUDE_DIRS} )

# bzip2
message (STATUS "Checking for bzip2")
find_package (BZip2 REQUIRED)
include_directories( ${BZIP_INCLUDE_DIRS} )

# boost
message (STATUS "Checking for Boost")
find_package(Boost 1.57.0 REQUIRED)
include_directories( ${Boost_INCLUDE_DIRS} )

# snappy
ExternalProject_Add(
  project_snappy
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/snappy
  BUILD_IN_SOURCE 1
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/snappy
  CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/snappy/autogen.sh COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/snappy/configure
  BUILD_COMMAND $(MAKE)
  INSTALL_COMMAND ""
)

ExternalProject_Get_Property(project_snappy install_dir)
add_library(snappy STATIC IMPORTED)
include_directories(${install_dir})

set(snappy_location ${install_dir}/.libs/libsnappy.a)
message("snappy target location is ${snappy_location}")
set_property(TARGET snappy PROPERTY IMPORTED_LOCATION ${snappy_location})
add_dependencies(snappy project_snappy)

## rocksdb
# Needed so rocksdb finds snappy
set(EXTRA_LDFLAGS "-L${CMAKE_CURRENT_SOURCE_DIR}/snappy/.libs/")
set(EXTRA_CFLAGS "-I${CMAKE_CURRENT_SOURCE_DIR}/snappy/")
set(ENV{EXTRA_LDFLAGS} ${EXTRA_LDFLAGS})
set(ENV{EXTRA_CFLAGS} ${EXTRA_CFLAGS})
set(ENV{CFLAGS} ${EXTRA_CFLAGS})

#message("Extra cflags = ${EXTRA_CFLAGS}")
#message("Extra ldflags = ${EXTRA_LDFLAGS}")

ExternalProject_Add(
  project_rocksdb
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/rocksdb
  BUILD_IN_SOURCE 1
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/rocksdb
  CONFIGURE_COMMAND ""
  BUILD_COMMAND PORTABLE=1;CFLAGS=${EXTRA_CFLAGS} $(MAKE) -C ${CMAKE_CURRENT_SOURCE_DIR}/rocksdb static_lib
  INSTALL_COMMAND ""
)

ExternalProject_Get_Property(project_rocksdb install_dir)
add_library(rocksdb STATIC IMPORTED)
include_directories(${install_dir}/include)
set_property(TARGET rocksdb PROPERTY IMPORTED_LOCATION ${install_dir}/librocksdb.a)
add_dependencies(project_rocksdb project_snappy snappy)
add_dependencies(rocksdb project_rocksdb project_snappy snappy)

# htslib
ExternalProject_Add(
  project_htslib
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/htslib
  BUILD_IN_SOURCE 1
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/htslib
  CONFIGURE_COMMAND ""
  BUILD_COMMAND $(MAKE) -C ${CMAKE_CURRENT_SOURCE_DIR}/htslib
  INSTALL_COMMAND ""
)

ExternalProject_Get_Property(project_htslib install_dir)
add_library(htslib STATIC IMPORTED)
include_directories(${install_dir})
include_directories(${install_dir}/htslib)

set(htslib_location ${install_dir}/libhts.a)
message("htslib target location is ${htslib_location}")
set_property(TARGET htslib PROPERTY IMPORTED_LOCATION ${htslib_location})
add_dependencies(htslib project_htslib)

# ExternalProject_Get_Property(htslib source_dir)
# set(htslib_INCLUDE_DIRS ${source_dir} CACHE INTERNAL "Path to include folder for htslib.")
# set(HTSlib_LIBRARY ${source_dir}/libhts.a CACHE INTERNAL "Path to htslib library files.")

# CMake SeqAn
# add_subdirectory(seqan)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/seqan/include)

# message("SeqAn include path is: ${SEQAN_INCLUDE_PATH}")
# include_directories(${SEQAN_INCLUDE_DIRS})
# include_directories(${SEQAN_INCLUDE_PATH})
# set(LIBS ${LIBS} ${SEQAN_LIBRARIES})

# Executables go to a "bin" folder.
set(EXECUTABLE_OUTPUT_PATH ../bin)

# Include public header files
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

## Compile gyper source files
message(STATUS "Compiling gyper's source files")
add_subdirectory(src)

# message("Including HTSlib: " ${htslib_INCLUDE_DIRS})
# include_directories(${htslib_INCLUDE_DIRS})

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

