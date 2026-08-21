cmake_minimum_required(VERSION 3.20)
project(Hemo1D LANGUAGES CXX)
 
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
 
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
 
option(HEMO1D_BUILD_TESTS "Build Hemo1D tests" ON)
 
include(FetchContent)
 
FetchContent_Declare(
    nlohmann_json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    GIT_TAG v3.11.3
)
FetchContent_MakeAvailable(nlohmann_json)
 
find_package(OpenMP QUIET)
 
add_subdirectory(src)
add_subdirectory(apps)
 
if(HEMO1D_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()