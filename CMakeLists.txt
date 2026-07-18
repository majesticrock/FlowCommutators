# FlowCommutators/CMakeLists.txt

cmake_minimum_required(VERSION 3.30)

project(FlowCommutators VERSION 1.0)

# Set local user-specific settings
# By default, this adds ~/usr/local/ and ../../.mrock/ to the cmake and include paths
# If additional configurations/includes are needed, provide them in cmake/user_settings.cmake
# Extra include dirs are stored in the cmake variable EXTRA_INCLUDE_DIRS
include(cmake/include_settings.cmake)

# Build type settings; defaults to Release
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

#
# Library specifics
#

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
