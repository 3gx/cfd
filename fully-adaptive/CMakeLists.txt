# CMake project
cmake_minimum_required(VERSION 2.8)
project("project_sigma" C CXX)

# Options
option(ENABLE_SAMPLES "Enable compilation of sample programs" OFF)
option(ENABLE_TESTS "Enable compilation of the Google tests" ON)



# Compiler-version check
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9)
    message(FATAL_ERROR "GCC version must be at least 4.9 (for C++14)")
  endif()
endif()

# C++14 compiler settings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++1y")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++1y -Wall -mavx2 -march=core-avx2")

# Package scripts location
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# The includes
include_directories(
  ${project_sigma_SOURCE_DIR}/include
  )

if (ENABLE_SAMPLES)
endif()

if (ENABLE_TESTS)
  add_subdirectory(external/gtest) 
  enable_testing()
  include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

  add_executable(unit_tests 
    unit_tests/m2n.cpp
    )

  add_executable(m2n
    src/m2n.cpp
    )

  target_link_libraries(unit_tests gtest gtest_main)
  target_link_libraries(m2n lapack)

  # Adds the tests
  add_test(
    name unit_tests
    command unit_tests
    )




endif()

