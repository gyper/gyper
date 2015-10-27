cmake_minimum_required(VERSION 2.8)

# Create a Gyper library
set(gyper_SOURCE_FILES
   graph.cpp
   graph_align.cpp
   graph_builder.cpp
   graph_io.cpp
   graph_kmerify.cpp
)
# file(GLOB gyper_SOURCE_FILES *.cpp)
# list(REMOVE_ITEM gyper_SOURCE_FILES main.cpp)
add_library(gyper_core ${gyper_SOURCE_FILES} )
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(gyper main.cpp)
find_package( Threads REQUIRED )
target_link_libraries(gyper gyper_core pthread ${ZLIB_LIBRARIES} ${Boost_LIBRARIES})

# Objectify Gyper
set(gyper_oo_SOURCE_FILES
   partial_order_graph.cpp
)

add_library(gyper_oo_core ${gyper_oo_SOURCE_FILES} )
target_link_libraries(gyper gyper_core pthread ${ZLIB_LIBRARIES} ${Boost_LIBRARIES})

# Add the executable
target_link_libraries(gyper ${LIBS})

# Installs gyper to /usr/local/bin using [sudo] make install
install (TARGETS gyper DESTINATION bin)




