cmake_minimum_required(VERSION 2.8)

# Create a Gyper library
set(gyper_SOURCE_FILES
   graph.cpp
   graph_align.cpp
   graph_builder.cpp
   graph_io.cpp
   graph_kmerify.cpp
)

# message (STATUS "Checking for snappy")
# ExternalProject_Get_Property(snappy install_dir)
# set(snappy_ROOT ${install_dir})
# 
# find_package(snappy REQUIRED)
# message("snappy_INCLUDE_DIRS = ${snappy_INCLUDE_DIRS}")
# message("snappy_LIBRARIES = ${snappy_LIBRARIES}")
# include_directories(${snappy_INCLUDE_DIRS})
# link_directories(${snappy_LIBRARIES})
# 
# # Make rocksdb
# set(ENV{EXTRA_CFLAGS} "-I${snappy_INCLUDE_DIRS}")
# set(ENV{EXTRA_CXXFLAGS} "-I${snappy_INCLUDE_DIRS}")
# set(ENV{EXTRA_LDFLAGS} "-L${snappy_LIBRARIES}")
# 
# ExternalProject_Get_Property(rocksdb install_dir)
# include_directories(${install_dir}/rocksdb/include)
# 
# set(ROCKSDB_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/rocksdb")
# find_library(ROCKSDB_LIBRARY NAMES rocksdb PATHS ${ROCKSDB_ROOT})
# ExternalProject_Get_Property(rocksdb source_dir)
# set(ROCKSDB_INCLUDE_DIR ${source_dir}/include CACHE INTERNAL "Path to include folder for rocksdb.")
# set(ROCKSDB_LIBRARIES ${source_dir}/librocksdb.a CACHE INTERNAL "Path to rocksdb library files.")
# include_directories(${ROCKSDB_INCLUDE_DIR})
# link_directories(${ROCKSDB_LIBRARIES})
# 
# message("ROCKSDB_ROOT = ${ROCKSDB_ROOT}")
# message("ROCKSDB_INCLUDE_DIR = ${ROCKSDB_INCLUDE_DIR}")
# message("ROCKSDB_LIBRARIES = ${ROCKSDB_LIBRARIES}")
# 
add_library(gyper_core ${gyper_SOURCE_FILES} )
add_dependencies(gyper_core project_rocksdb project_snappy project_htslib seqan)

add_executable(gyper main.cpp)
add_dependencies(gyper project_rocksdb project_snappy project_htslib htslib rocksdb snappy seqan)
find_package( Threads REQUIRED )
target_link_libraries(gyper gyper_core htslib rocksdb snappy pthread ${ZLIB_LIBRARIES} ${Boost_LIBRARIES})

# Objectify Gyper
set(gyper_oo_SOURCE_FILES
   constructor.cpp
   options.cpp
   graph_creation.cpp
)

add_library(gyper_oo_core ${gyper_oo_SOURCE_FILES} )
add_dependencies(gyper_oo_core project_rocksdb project_snappy project_htslib seqan)
target_link_libraries(gyper gyper_core gyper_oo_core htslib pthread ${ZLIB_LIBRARIES} ${Boost_LIBRARIES})

# Add the executable
target_link_libraries(gyper ${LIBS})

# Installs gyper to /usr/local/bin using [sudo] make install
install (TARGETS gyper DESTINATION bin)
