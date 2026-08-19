set(HEMO1D_TEST_SOURCES
    network_test.cpp
    test_linear_algebra.cpp
)

if(NOT HEMO1D_TEST_SOURCES)
    message(STATUS "hemo1d: no test sources yet. Skipping hemo1d_tests")
    return()
endif()

include(FetchContent)
 
FetchContent_Declare(
    catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.5.4
)
FetchContent_MakeAvailable(catch2)
 
add_executable(hemo1d_tests ${HEMO1D_TEST_SOURCES})

target_link_libraries(hemo1d_tests PRIVATE hemo1d_core Catch2::Catch2WithMain)
 
target_compile_definitions(hemo1d_tests PRIVATE
    HEMO1D_EXAMPLES_DIR="${PROJECT_SOURCE_DIR}/examples/networks"
)
 
list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(Catch)
catch_discover_tests(hemo1d_tests)