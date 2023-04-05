#include "galerkin_projection_rom.hpp"
#include "burgers_rewienski.hpp"
#include <iostream>
#include "misc/eigen_utils.hpp"

GalerkinProjectionROM::GalerkinProjectionROM(
    Snapshots &snapshot_matrix,
    const double &reduced_iteration_tolerance,
    const unsigned int &max_reduced_iterations)
    : snapshot_matrix(snapshot_matrix),
    reduced_iteration_tolerance(reduced_iteration_tolerance),
    max_reduced_iterations(max_reduced_iterations)
{
    basis = snapshot_matrix.get_basis();
    reference_state = snapshot_matrix.get_reference_state();
}

void GalerkinProjectionROM::steady_state(const float &b)
{
    std::cout << "\n-------------------\nSolving ROM for b = " << b << "\n-------------------" << std::endl;

    rom_solution.resize(snapshot_matrix.fom_solver.nx);
    for (int i = 0; i < snapshot_matrix.fom_solver.nx; i++)  // Setting intial condition
    {
        rom_solution(i) = snapshot_matrix.fom_solver.bc_spec.second[0] ;
    }
    rom_solution = reference_state + basis * (basis.transpose() * (rom_solution - reference_state));  // Projecting initial guess

    Eigen::VectorXd residual = Eigen::VectorXd::Zero(snapshot_matrix.fom_solver.nx);  // Residual in each cell for the current iteration
    Eigen::MatrixXd jacobian = Eigen::MatrixXd::Ones(snapshot_matrix.fom_solver.nx, snapshot_matrix.fom_solver.nx);

    double reduced_residual = 1;
    unsigned int itr = 0;
    while (reduced_residual > reduced_iteration_tolerance && itr < max_reduced_iterations)
    {
        std::cout << "\nReduced Iteration: " << itr
                  << "\tCurrent Reduced Residual: " << reduced_residual
                  << std::endl;

        snapshot_matrix.fom_solver.evaluate_residual(rom_solution, residual, b);  // Full order residual using approximate solution
        snapshot_matrix.fom_solver.evaluate_jacobian(rom_solution, jacobian, b);  // Full order Jacobian using approximate solution
        reduced_residual = reduced_residual_L2_norm(jacobian, residual);

        // Computing the search direction
        Eigen::MatrixXd reduced_LHS = basis.transpose() * jacobian.transpose() * jacobian * basis;
        Eigen::VectorXd reduced_RHS = -1 * basis.transpose() * jacobian.transpose() * residual;
        Eigen::VectorXd search_direction = reduced_LHS.householderQr().solve(reduced_RHS);

        // Computing step size using backtracking algorithm
        double step_length = line_search(rom_solution, reduced_residual, search_direction, b);

        rom_solution = rom_solution + basis * (step_length * search_direction);

        itr++;
    }
}

double GalerkinProjectionROM::reduced_residual_L2_norm(const Eigen::MatrixXd &jacobian, const Eigen::VectorXd &residual)
{
    Eigen::VectorXd reduced_residual = basis.transpose() * jacobian.transpose() * residual;
    return reduced_residual.norm();
}

double GalerkinProjectionROM::line_search(
    const Eigen::VectorXd &old_rom_solution,
    const double &old_reduced_residual,
    const Eigen::VectorXd &search_direction,
    const float &b)
{
    double step_length = 100;
    float step_length_reduction_factor = 0.95;
    float min_step_length = 1e-5;
    unsigned int max_line_seach_itr = 100;
    float residual_reduction_tolerance = 1;
    bool switch_search_direction = false;

    Eigen::VectorXd new_residual = Eigen::VectorXd::Zero(snapshot_matrix.fom_solver.nx);
    Eigen::MatrixXd new_jacobian = Eigen::MatrixXd::Ones(snapshot_matrix.fom_solver.nx, snapshot_matrix.fom_solver.nx);
    Eigen::VectorXd new_rom_solution;

    double new_reduced_residual = 1e6;
    unsigned int itr = 0;

    std::cout << "Begining line search" << std::endl;

    while (new_reduced_residual > old_reduced_residual * residual_reduction_tolerance
           && itr < max_line_seach_itr && std::abs(step_length) > min_step_length)
    {
        new_rom_solution = old_rom_solution + basis * (step_length * search_direction);
        snapshot_matrix.fom_solver.evaluate_residual(new_rom_solution, new_residual, b);
        snapshot_matrix.fom_solver.evaluate_jacobian(new_rom_solution, new_jacobian, b);
        new_reduced_residual = reduced_residual_L2_norm(new_jacobian, new_residual);
        itr++;

        std::cout << "\tForwards line search iteration: " << itr
                  << "\tStep length: " << step_length
                  << "\tNew reduced residual: " << new_reduced_residual
                  << std::endl;

        step_length *= step_length_reduction_factor;
    }
    if (itr == max_line_seach_itr && new_reduced_residual > old_reduced_residual) 
    {
        switch_search_direction = true;
        std::cout << "\n\tReached max line search iterations, reduced residual decrease not sufficient. Switching search direction\n" << std::endl;
    }

    if (switch_search_direction)
    {
        unsigned int itr = 0;
        step_length = -100;
        new_reduced_residual = 1e6;
        while (new_reduced_residual > old_reduced_residual * residual_reduction_tolerance
               && itr < max_line_seach_itr && std::abs(step_length) > min_step_length)
        {
            new_rom_solution = rom_solution + basis * (step_length * search_direction);
            snapshot_matrix.fom_solver.evaluate_residual(new_rom_solution, new_residual, b);
            new_reduced_residual = new_residual.norm();
            itr++;

            std::cout << "\tBackwards line search iteration: " << itr
                      << "\tStep length: " << step_length
                      << "\tNew reduced residual: " << new_reduced_residual
                      << std::endl;

            step_length *= step_length_reduction_factor;
        }
        if (itr == max_line_seach_itr && new_reduced_residual > old_reduced_residual) 
        {
            std::cout << "Line search failed, reduced residual did not decrease. Setting step length to 0.\n" << std::endl;
            return 0;
        }
    }
    return step_length;
}

