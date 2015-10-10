cmake_minimum_required(VERSION 2.8.8)

project(test_gyper CXX)

message("Compiling testing suite")

file(GLOB gyper_TEST_FILES *.cpp)

add_executable(test_gyper ${gyper_TEST_FILES})
target_link_libraries(test_gyper gyper_lib pthread ${ZLIB_LIBRARIES} ${Boost_LIBRARIES})
