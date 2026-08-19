set(HEMO1D_CORE_SOURCES
    core/network.cpp
    core/node.cpp
    core/vessel.cpp
    core/linear_algebra.cpp
)

if(NOT HEMO1D_CORE_SOURCES)
    message(STATUS "hemo1d: no library sources yet. Skipping hemo1d_core")
    return()
endif()

add_library(hemo1d_core ${HEMO1D_CORE_SOURCES})

target_include_directories(hemo1d_core PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)

target_compile_features(hemo1d_core PUBLIC cxx_std_20)

target_link_libraries(hemo1d_core PUBLIC nlohmann_json::nlohmann_json)

if(OpenMP_CXX_FOUND)
    target_link_libraries(hemo1d_core PUBLIC OpenMP::OpenMP_CXX)
endif()