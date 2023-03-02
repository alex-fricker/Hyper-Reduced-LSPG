add_subdirectory("tests")

find_package(Eigen3 3.3 REQUIRED)

add_executable(main main.cpp)

add_library(burgers_rewienski burgers_rewienski.cpp)

target_link_libraries(
    main 
    tests
    burgers_rewienski
    Eigen3::Eigen
)

