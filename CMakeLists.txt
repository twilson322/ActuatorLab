cmake_minimum_required(VERSION 3.22)
project(actuator_lab VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(ACTLAB_SANITIZE "Build with AddressSanitizer + UndefinedBehaviorSanitizer" OFF)
option(ACTLAB_TSAN "Build with ThreadSanitizer" OFF)
option(ACTLAB_BUILD_TESTS "Build unit tests" ON)

add_compile_options(-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wno-sign-conversion)
if(ACTLAB_SANITIZE)
  add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer -g)
  add_link_options(-fsanitize=address,undefined)
endif()
if(ACTLAB_TSAN)
  add_compile_options(-fsanitize=thread -g -O1)
  add_link_options(-fsanitize=thread)
endif()

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(Threads REQUIRED)

add_library(actlab
  src/plant.cpp
  src/controllers.cpp
  src/lqg.cpp
  src/drive.cpp
  src/experiments.cpp)
target_include_directories(actlab PUBLIC include)
target_link_libraries(actlab PUBLIC Eigen3::Eigen Threads::Threads)

add_executable(actlab_sim apps/actlab_sim.cpp)
target_link_libraries(actlab_sim PRIVATE actlab)

add_executable(rt_loop_bench apps/rt_loop_bench.cpp)
target_link_libraries(rt_loop_bench PRIVATE actlab)

if(ACTLAB_BUILD_TESTS)
  enable_testing()
  find_package(GTest REQUIRED)
  add_executable(actlab_tests
    tests/test_foc.cpp
    tests/test_filters.cpp
    tests/test_rls.cpp
    tests/test_plant.cpp
    tests/test_spsc.cpp
    tests/test_closed_loop.cpp)
  target_link_libraries(actlab_tests PRIVATE actlab GTest::gtest_main)
  include(GoogleTest)
  gtest_discover_tests(actlab_tests)
endif()
