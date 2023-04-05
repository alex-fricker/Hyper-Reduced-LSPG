#include "../burgers_rewienski.hpp"
#include <iostream>
#include "../misc/eigen_utils.hpp"


int run_1D_burgers_rewienski_solution()
{
    const int nx=200;
    const float x0=0;
    const float x1=75;

    std::vector<float> bc_vals(1, 1);
    std::vector<float> ic_vals(1, 1);
    std::pair<std::string, std::vector<float>> BC("Constant", bc_vals);
    std::pair<std::string, std::vector<float>> IC("Constant", ic_vals);

    BurgersRewienski solver = BurgersRewienski(nx, x0, x1, BC, IC);
    std::vector<float> b_range {0.1};

    for (float b:b_range)
    {
        std::cout << "\n===========================\nRunning for b=" << b << "\n===========================" << std::endl;

        Eigen::VectorXd soln = solver.solve(b);
        Eigen::VectorXd residual = std2eigen_vector(solver.get_residual_history());
        write_vector("burgers_rewienski_b_" + std::to_string(b) + ".txt", soln);
        write_vector("burgers_rewienski_residual_b_" + std::to_string(b) + ".txt", residual);
    }
    return 0;
}
