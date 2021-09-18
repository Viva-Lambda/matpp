cmake_minimum_required (VERSION 3.0.2)
project("matpp")

# include ctest
include(CTest)

#
enable_testing()

set (CMAKE_BUILD_TYPE "Debug")

set (CMAKE_CXX_FLAGS "-std=c++14")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -ggdb -Wall -Wextra -ldl")
set(CMAKE_EXPORT_COMPILE_COMMANDS true)


# include test suite
include_directories("./include/")

# test dir
set (TEST_DIR "${PROJECT_SOURCE_DIR}/tests")

set (TEST_MAIN "${PROJECT_SOURCE_DIR}/test.c")

file(GLOB test_files "${TEST_DIR}/*.cpp")

foreach(test_file ${test_files})
    message(STATUS "test file ${test_file}")
    string(REPLACE "${TEST_DIR}/" "" file_no_parent "${test_file}")
    string(REPLACE ".cpp" ".out" file_exec_name "${file_no_parent}")
    add_executable(
        ${file_exec_name}  # executable name
        ${test_file}  # test file
        ${TEST_MAIN}  # test main entry point
    )
    message(STATUS "test file executable ${file_exec_name}")
# target_link_libraries(${file_exec_name} ${PROJECT_LIBS})
endforeach()

