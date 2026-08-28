# Let only one source to compile !!!!!!
# network_info network_info.cpp
add_executable(
    hemo1d_bench hemo1d_bench.cpp
)
target_link_libraries(hemo1d_bench PRIVATE hemo1d_core)