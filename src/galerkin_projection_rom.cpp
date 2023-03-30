#include "galerkin_projection_rom.hpp"
#include "burgers_rewienski.hpp"
#include <iostream>

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
    rom_solution.resize(snapshot_matrix.fom_solver.nx);
    for (int i = 0; i < snapshot_matrix.fom_solver.nx; i++)  // Setting intial condition
    {
        rom_solution(i) = snapshot_matrix.fom_solver.bc_spec.second[0] ;
    }
    rom_solution = reference_state + basis * (basis.transpose() * (rom_solution - reference_state));  // Projecting initial guess

    double reduced_residual = 1;
    unsigned int itr = 0;
    while (reduced_residual > reduced_iteration_tolerance && itr < max_reduced_iterations)
    {
        Eigen::VectorXd residual = snapshot_matrix.fom_solver.rom_residual(rom_solution, b);  // Full order residual using approximate solution
        Eigen::MatrixXd jacobian = snapshot_matrix.fom_solver.rom_jacobian(rom_solution, b);  // Full order Jacobian using approximate solution
        reduced_residual = reduced_residual_L2_norm(jacobian, residual);

        // Computing the search direction
        Eigen::MatrixXd reduced_LHS = basis.transpose() * jacobian.transpose() * jacobian * basis;
        Eigen::VectorXd reduced_RHS = -1 * basis.transpose() * jacobian.transpose() * residual;
        Eigen::VectorXd search_direction = reduced_LHS.householderQr().solve(reduced_RHS);

        // Computing step size using backtracking algorithm
        double step_length = line_search(rom_solution, search_direction, b);

        rom_solution += rom_solution + basis * (step_length * search_direction);

        itr++;
    }


}

double GalerkinProjectionROM::reduced_residual_L2_norm(const Eigen::MatrixXd &jacobian, const Eigen::VectorXd &residual)
{
    Eigen::MatrixXd reduced_residual = basis.transpose() * jacobian.transpose() * residual;
    return reduced_residual.norm();
}

double GalerkinProjectionROM::line_search(
    const Eigen::VectorXd &rom_solution,
    const Eigen::VectorXd &search_direction,
    const float &b)
{
    double step_length = 1;
    float step_length_reduction_factor = 0.5;
    unsigned int max_line_seach_itr = 10;
    float residual_reduction_tolerance = 1;
    bool switch_search_direction = false;


    double old_reduced_residual = snapshot_matrix.fom_solver.rom_residual(rom_solution, b).norm();
    double new_reduced_residual = 1;
    unsigned int itr = 0;
    while (new_reduced_residual > old_reduced_residual * residual_reduction_tolerance)
    {
        step_length *= step_length_reduction_factor;
        Eigen::VectorXd new_rom_solution = rom_solution + step_length * search_direction;
        new_reduced_residual = snapshot_matrix.fom_solver.rom_residual(new_rom_solution, b).norm();

        if (itr == max_line_seach_itr && new_reduced_residual > old_reduced_residual) { switch_search_direction = true; }

    }

    if (switch_search_direction)
    {
        unsigned int itr = 0;
        step_length = -1;
        while (new_reduced_residual > old_reduced_residual * residual_reduction_tolerance)
        {
            step_length *= step_length_reduction_factor;
            new_reduced_residual = snapshot_matrix.fom_solver.rom_residual(rom_solution + step_length * search_direction, b).norm();

            if (itr == max_line_seach_itr && new_reduced_residual > old_reduced_residual) 
            {
                std::cout << "Line search failed, reduced residual did not decrease. Setting step length to 0.\n" << std::endl;
                return 0;
            }
        }
    }
    return step_length;
}

