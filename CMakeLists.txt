# FlowCommutators/CMakeLists.txt

cmake_minimum_required(VERSION 3.30)

project(FlowCommutators VERSION 1.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_PREFIX_PATH "$ENV{HOME}/usr/local")

add_executable(FlowCommutators 
    sources/FlowCommutators.cpp
)

target_compile_features(FlowCommutators
    INTERFACE
        cxx_std_20
)

find_package(mrock 
    REQUIRED 
    COMPONENTS 
        symbolic_operators 
)
target_link_libraries(FlowCommutators PUBLIC ${mrock_LIBRARIES})
