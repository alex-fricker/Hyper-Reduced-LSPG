add_subdirectory("tests")

find_package(Eigen3 3.3 REQUIRED)

add_executable(main main.cpp)

add_library(burgers_rewienski burgers_rewienski.cpp)
add_library(snapshots snapshots.cpp)
add_library(galerkin_projection_rom galerkin_projection_rom.cpp)

target_link_libraries(
    main 
    tests
    burgers_rewienski
    Eigen3::Eigen
    snapshots
    galerkin_projection_rom
)

