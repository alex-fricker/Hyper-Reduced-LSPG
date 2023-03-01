project(
HRLSPG
DESCRIPTION "Hyper-reduced least squares Petrov Galerkin projection for the Burgers Rewienski equation."
LANGUAGES CXX)

CMAKE_MINIMUM_REQUIRED(VERSION 3.3.0)
set(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
    
set(WARNING_CXX_FLAGS "-Wall -Werror -pedantic -pedantic-errors")
set(WARNING_CXX_FLAGS "${WARNING_CXX_FLAGS} -Wextra -Wunreachable-code")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${WARNING_CXX_FLAGS} -march=native -std=gnu++17")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

add_subdirectory(src)


