#ifndef __GALERKIN_PROJECTION_ROM__
#define __GALERKIN_PROJECTION_ROM__

#include "../Libraries/eigen/Eigen/Dense"
// #include "../Libraries/CoDiPack/include/codi.hpp"
#include "snapshots.hpp"

class GalerkinProjectionROM
{
public:
    // Constructor
    GalerkinProjectionROM(
        Snapshots &snapshot_matrix,
        const double &reduced_iteration_tolerance,
        const unsigned int &max_reduced_iterations);

    // Destructor
    ~GalerkinProjectionROM() {};

    // Run the ROM to solve the approximate solution
    void steady_state(const float &b);

    Eigen::VectorXd rom_solution;  // Approximate solution at iteration k
    
protected:
    Snapshots snapshot_matrix;
    Eigen::MatrixXd basis;
    Eigen::VectorXd reference_state;

    // Compute the L2 norm of the reduced residual
    double reduced_residual_L2_norm(const Eigen::MatrixXd &jacobian, const Eigen::VectorXd &residual);

    double line_search(
        const Eigen::VectorXd &rom_solution,
        const double &old_reduced_residual,
        const Eigen::VectorXd &search_direction,
        const float &b);

    const double reduced_iteration_tolerance;
    const double max_reduced_iterations;
};

#endif
