cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}")

add_subdirectory(rvlm-common)
add_subdirectory(netcdf-cxx4)
add_subdirectory(EasyCL)

add_executable(fdtd
   src/main.cc
   #src/read_grid.cc
   src/yee_grid.hh
   src/resistive_source.hh
   src/kernels.cl
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

target_include_directories(fdtd PRIVATE "${NETCDF_INCLUDE_DIRS}" "EasyCL")
target_link_libraries(fdtd rvlm-common netcdf-cxx4 EasyCL)
