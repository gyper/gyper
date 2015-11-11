cmake_minimum_required(VERSION 2.8)

project(test_gyper CXX)


set(gyper_TEST_FILES test.cpp test_gyper_constructor.cpp test_seqanhts.cpp)

add_executable(test_gyper ${gyper_TEST_FILES})
add_dependencies(test_gyper project_rocksdb project_snappy project_htslib seqan)

# gyper_core
target_link_libraries(test_gyper gyper_oo_core htslib rocksdb snappy pthread ${ZLIB_LIBRARIES} ${BZIP2_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_rocksdb test_rocksdb.cpp)
target_link_libraries(test_rocksdb ${ZLIB_LIBRARIES} ${BZIP2_LIBRARIES} rocksdb snappy)

add_test(NAME RunTests COMMAND test_gyper)

add_test(NAME ListTests COMMAND test_gyper --list-tests)
set_tests_properties(ListTests PROPERTIES PASS_REGULAR_EXPRESSION "[0-9]+ test cases")

add_test(NAME ListTags COMMAND test_gyper --list-tags)
set_tests_properties(ListTags PROPERTIES PASS_REGULAR_EXPRESSION "[0-9]+ tags")
