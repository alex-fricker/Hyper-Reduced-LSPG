#include "snapshots_test.hpp"
#include "../snapshots.hpp"
#include "../burgers_rewienski.hpp"
#include "../misc/eigen_utils.cpp"
#include <string>

int run_snapshots_test()
{
    // Snapshots parameters
    const int number_snapshots = 3;
    const std::pair<double, double> b_range (0.01, 0.1);
    const double evaluation_time = 1;

    // Build FOM solver
    const int nx=1024;
    const float x0=0;
    const float x1=100;
    std::vector<float> bc_vals(1, 1);
    std::vector<float> ic_vals(1, 1);
    std::pair<std::string, std::vector<float>> BC("Constant", bc_vals);
    std::pair<std::string, std::vector<float>> IC("Constant", ic_vals);
    BurgersRewienski fom = BurgersRewienski(nx, x0, x1, BC, IC);

    // Build snapshot matrix
    Snapshots snapshot_matrix = Snapshots(fom);
    snapshot_matrix.build_snapshot_matrix(number_snapshots, b_range, evaluation_time);
    snapshot_matrix.compute_pod_basis();

    const std::string basis_name = "basis_matrix_" + std::to_string(number_snapshots) + "_basis_vectors";
    const std::string snap_name = "snapshot_matrix_" + std::to_string(number_snapshots) + "_snapshots";
    const std::string points_name = "snapshot_points_" + std::to_string(number_snapshots) + "_snapshots";
    const std::string residual_name = "snapshot_residuals_" + std::to_string(number_snapshots) + "_snapshots";
    write_matrix(snap_name, snapshot_matrix.get_snapshot_matrix());
    write_matrix(basis_name, snapshot_matrix.get_basis());
    write_vector(points_name, snapshot_matrix.get_snapshot_points());
    write_vector(residual_name, snapshot_matrix.get_snapshot_residuals());

    return 0;
}
