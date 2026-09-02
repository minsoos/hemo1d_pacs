include(FetchContent)
FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG v2.13.6
)
FetchContent_MakeAvailable(pybind11)
 
pybind11_add_module(hemo1d_py bindings.cpp)
target_link_libraries(hemo1d_py PRIVATE hemo1d_core hemo1d_couplings)
 
set_target_properties(hemo1d_py PROPERTIES OUTPUT_NAME "hemo1d")
 
install(TARGETS hemo1d_py LIBRARY DESTINATION . RUNTIME DESTINATION .)