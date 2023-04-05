#include "../galerkin_projection_rom.hpp"
#include "galerkin_projection_test.hpp"
#include "../../Libraries/eigen/Eigen/Dense"

int run_galerkin_projection_test()
{
    // Snapshots parameters
    const int number_snapshots = 100;
    const std::pair<double, double> b_range (0.01, 0.1);

    // Build FOM solver
    const int nx=200;
    const float x0=0;
    const float x1=75;
    std::vector<float> bc_vals(1, 1);
    std::vector<float> ic_vals(1, 1);
    std::pair<std::string, std::vector<float>> BC("Constant", bc_vals);
    std::pair<std::string, std::vector<float>> IC("Constant", ic_vals);
    BurgersRewienski fom = BurgersRewienski(nx, x0, x1, BC, IC);

    // Build snapshot matrix
    Snapshots snapshot_matrix = Snapshots(fom);
    snapshot_matrix.build_snapshot_matrix(number_snapshots, b_range);
    snapshot_matrix.compute_pod_basis();

    // Build LSPG
    double reduced_iteration_tolerance = 1e-1;
    unsigned int max_reduced_iterations = 1000;
    GalerkinProjectionROM rom = GalerkinProjectionROM(snapshot_matrix, reduced_iteration_tolerance, max_reduced_iterations);
    float rom_point = 0.055;
    rom.steady_state(rom_point);
    Eigen::VectorXd rom_solution = rom.rom_solution;


    return 0;
}
