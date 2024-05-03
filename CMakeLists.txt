cmake_minimum_required(VERSION 3.5)
project(random-orbits)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(PkgConfig REQUIRED)

pkg_search_module(flint REQUIRED flint>=3.0)
include_directories(${flint_INCLUDE_DIRS})
link_libraries(${flint_LDFLAGS})

add_executable(random-cubic random-cubic.cpp)
add_executable(random-quadratic random-quadratic.cpp)
