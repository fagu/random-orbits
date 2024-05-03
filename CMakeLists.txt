cmake_minimum_required(VERSION 3.5)
project(random-orbits)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(PkgConfig REQUIRED)

pkg_search_module(flint flint>=3.0)
if(${flint_FOUND})
	message("Pkgconf found flint module.")
	include_directories(${flint_INCLUDE_DIRS})
	link_libraries(${flint_LDFLAGS})
else()
	message("Pkgconf did not find flint module.")
	
	find_library(flint_old flint REQUIRED)
	link_libraries(${flint_old})
	
	find_library(gmp gmp REQUIRED)
	link_libraries(${gmp})
	
	find_library(arb_old flint-arb)
	if(${arb_old})
		link_libraries(${arb_old})
	else()
		message(WARNING "Could not find arb library.")
	endif()
endif()


add_executable(random-cubic random-cubic.cpp)
add_executable(random-quadratic random-quadratic.cpp)
