# FlowCommutators/CMakeLists.txt

cmake_minimum_required(VERSION 3.16)

project(FlowCommutators VERSION 1.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_PREFIX_PATH "~/usr/local")

add_executable(FlowCommutators 
    sources/FlowCommutators.cpp
)

# Check if using GCC and adjust C++ standard flag if necessary
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(FlowCommutators PRIVATE
        -Wall
        -Wno-sign-compare
        -march=native
        -O3
    )
else()
    message(FATAL_ERROR "Unsupported compiler ${CMAKE_CXX_COMPILER_ID}. Only GCC|Clang is supported.")
endif()

find_package(mrock 
    REQUIRED 
    COMPONENTS 
        symbolic_operators 
)
target_link_libraries(FlowCommutators PUBLIC ${mrock_LIBRARIES})
