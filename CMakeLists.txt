cmake_minimum_required(VERSION 2.8)
project(fdtd-experiments)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}")

add_subdirectory(rvlm-common)

add_executable(fdtd
   src/calc_field.hh
   src/resistive_source.hh
   src/CurrentSource.hh
   src/yee_grid.hh
   src/Dimensions.hh
   src/Constants.hh
   src/BoostUnitsHelpers.hh
   src/main.cc)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")

target_link_libraries(fdtd rvlm-common -lmpfr)
