set(HEMO1D_CORE_SOURCES
    physics/tube_law.cpp
    physics/characteristics.cpp
    physics/flux.cpp
    core/network.cpp
    core/node.cpp
    core/vessel.cpp
    io/time_series.cpp
    core/linear_algebra.cpp
    dg/quadrature.cpp
    dg/nodal_basis.cpp
    dg/reference_element.cpp
    dg/element.cpp
    dg/mesh.cpp
    io/network_parser.cpp
    physics/compatibility.cpp
    physics/spatial_operator.cpp
    physics/cfl.cpp
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