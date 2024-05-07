cmake_minimum_required(VERSION 3.18)
project(cylinder)

find_package(Open3D REQUIRED)

# Add the executable
add_executable(meshing-test)
target_sources(meshing-test PRIVATE main.cpp)

target_link_libraries(meshing-test PRIVATE Open3D::Open3D)
