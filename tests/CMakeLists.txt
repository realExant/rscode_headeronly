cmake_minimum_required(VERSION 3.21)

project(tests LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

Include(FetchContent)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.0.1 # or a later release
)

FetchContent_MakeAvailable(Catch2)

add_executable(tst_rscode tst_rscode.cpp)
target_link_libraries(tst_rscode PRIVATE Catch2::Catch2 Catch2::Catch2WithMain)

include(CTest)
include(Catch)
catch_discover_tests(tst_rscode)

