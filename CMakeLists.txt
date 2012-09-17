#===- CMakeLists.txt --------------------------------------------------------===
#
# Color
#
# This file is distributed under the Simplified BSD License. See LICENSE.TXT
# for details.
#
#===-------------------------------------------------------------------------===
#
# This is the make CMake file for Color
#
#===-------------------------------------------------------------------------===

cmake_minimum_required(VERSION 2.8)

project(Color)

set(Boost_USE_STATIC_LIBS TRUE)
set(Boost_ADDITIONAL_VERSIONS "1.50" "1.50.0" "1.44" "1.44.0" "1.42" "1.42.0" "1.41.0" "1.41" "1.40.0" "1.40" "1.39.0" "1.39" "1.38.0" "1.38" "1.37.0" "1.37" )
find_package(Boost COMPONENTS filesystem)
# Set up referencing of Boost
include_directories(${Boost_INCLUDE_DIR})
add_definitions(-D_SCL_SECURE_NO_WARNINGS)
# ${Boost_LIBRARIES}

add_executable(color
  color.cpp
  )

target_link_libraries(color ${Boost_LIBRARIES})
