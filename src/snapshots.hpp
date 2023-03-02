#ifndef __SNAPSHOTS__
#define __SNAPSHOTS__

#include "../Libraries/eigen/Eigen/Dense"
#include <vector>
#include "burgers_rewienski.hpp"


class Snapshots
{
public:
    // Constructor
    Snapshots(BurgersRewienski &fom_solver);

    // Destructor
    ~Snapshots() {};

    // Method to build the snapshot matrices
    void build_snapshot_matrix(
        const int n_snapshots,
        const std::pair<double, double> b_range,
        const double t_eval
    );

    // Method to compute the POD basis
    void compute_pod_basis();

    // Method for writing a Eigen matrix to a csv file
    void write_matrix();

    Eigen::MatrixXd get_snapshot_matrix();
    Eigen::MatrixXd get_basis();
    Eigen::VectorXd get_reference_state();

private:
    BurgersRewienski fom_solver;  // Solver for getting snapshot solutions
    Eigen::MatrixXd snapshot_matrix;  // Snapshot matrix
    Eigen::MatrixXd basis;  // POD basis
    Eigen::VectorXd reference_state;  // Reference state for snapshot matrix
};

#endif
