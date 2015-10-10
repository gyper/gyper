cmake_minimum_required(VERSION 2.8)

# Create a Gyper library
file(GLOB gyper_SOURCE_FILES *.cpp)
list(REMOVE_ITEM gyper_SOURCE_FILES main.cpp)
add_library(gyper_lib ${gyper_SOURCE_FILES} )
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(gyper main.cpp)
find_package( Threads REQUIRED )
target_link_libraries(gyper gyper_lib pthread ${ZLIB_LIBRARIES} ${Boost_LIBRARIES})


# Add the executable
target_link_libraries(gyper ${LIBS})

# Installs gyper to /usr/local/bin using [sudo] make install
install (TARGETS gyper DESTINATION bin)
