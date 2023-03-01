add_subdirectory("tests")

add_executable(main main.cpp)

add_library(burgers_rewienski burgers_rewienski.cpp)

target_link_libraries(
    main 
    tests
    burgers_rewienski
)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries (example Eigen3::Eigen)
