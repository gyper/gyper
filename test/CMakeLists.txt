cmake_minimum_required(VERSION 2.8)

project(test_gyper CXX)

message("Compiling testing suite")
message("Flags are: ${CMAKE_CXX_FLAGS}")

#file(GLOB gyper_TEST_FILES *.cpp)

set(gyper_TEST_FILES test.cpp test_graph_hla.cpp)

add_executable(test_gyper ${gyper_TEST_FILES})
# gyper_core
target_link_libraries(test_gyper gyper_oo_core hts pthread ${ZLIB_LIBRARIES} ${Boost_LIBRARIES})

add_test(NAME RunTests COMMAND test_gyper)

add_test(NAME ListTests COMMAND test_gyper --list-tests)
set_tests_properties(ListTests PROPERTIES PASS_REGULAR_EXPRESSION "[0-9]+ test cases")

add_test(NAME ListTags COMMAND test_gyper --list-tags)
set_tests_properties(ListTags PROPERTIES PASS_REGULAR_EXPRESSION "[0-9]+ tags")
