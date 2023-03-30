#include "tests/1D_burgers_rewienski_solution.hpp"
#include "tests/snapshots_test.hpp"
#include "tests/galerkin_projection_test.hpp"
#include<iostream>

int main()
{
    // int result = run_1D_burgers_rewienski_solution();
    // int result = run_snapshots_test();
    int result = run_galerkin_projection_test();

    return result;
}
