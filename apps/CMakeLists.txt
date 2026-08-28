add_executable(
    hemo1d_network_info network_info.cpp
    hemo1d_bench.cpp
)
target_link_libraries(hemo1d_network_info PRIVATE hemo1d_core)