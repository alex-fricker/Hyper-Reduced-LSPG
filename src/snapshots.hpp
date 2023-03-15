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

    const Eigen::MatrixXd get_snapshot_matrix() const; 
    const Eigen::MatrixXd get_basis() const;
    const Eigen::VectorXd get_reference_state() const;
    const Eigen::VectorXd get_snapshot_points() const;
    const Eigen::VectorXd get_snapshot_residuals() const;

private:
    BurgersRewienski fom_solver;  // Solver for getting snapshot solutions
    Eigen::MatrixXd snapshot_matrix;  // Snapshot matrix
    Eigen::MatrixXd basis;  // POD basis
    Eigen::VectorXd reference_state;  // Reference state for snapshot matrix
    Eigen::VectorXd snapshot_points;  // a values in first row, b values in second row
    Eigen::VectorXd snapshot_residuals;  // Vector to store the residual of each snapshot

    double halton_element(int i, int m);
    int i4vec_sum(int n, int a[]);
    int prime(int n);
};

#endif
