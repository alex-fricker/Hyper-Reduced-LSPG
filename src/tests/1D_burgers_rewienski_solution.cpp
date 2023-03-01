#include "../burgers_rewienski.hpp"
#include <iostream>


int run_1D_burgers_rewienski_solution()
{
    std::cout << "Inside 1D BR Sol" << std::endl;

    const int nx=1024;
    const float x0=0;
    const float x1=100;
    const float CFL=0.9;
    const float t1=0.5;

    std::vector<float> bc_vals(1, 1);
    std::vector<float> ic_vals(1, 1);
    std::pair<std::string, std::vector<float>> BC("Constant", bc_vals);
    std::pair<std::string, std::vector<float>> IC("Constant", ic_vals);

    BurgersRewienski solver = BurgersRewienski(nx, x0, x1, CFL, BC, IC);
    float b = 0.05;
    std::vector<float> soln = solver.solve(b, t1);
    solver.write_solution("burgers_rewienski_b_" + std::to_string(b), soln);

    std::vector<float> b_range {0.1, 0.05, 0.01};

    // for (float b:b_range)
    // {
    //     std::cout << "Running for b=" << b << std::endl;

    //     std::vector<float> soln = solver.solve(b, t1);
    //     solver.write_solution("burgers_rewienski_b_" + std::to_string(b), soln);
    // }
    return 0;
}
