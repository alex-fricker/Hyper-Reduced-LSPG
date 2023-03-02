#include <iostream>
#include <fstream>
#include "snapshots.hpp"
#include "../Libraries/eigen/Eigen/SVD"
#include "halton.hpp"


Snapshots::Snapshots(BurgersRewienski &fom_solver)
    : fom_solver(fom_solver)
{

}

void Snapshots::build_snapshot_matrix(
    const int n_snapshots, 
    const std::pair<double, double> b_range,
    const double t_eval
)
{
    // Determine points to evaluate snapshot matrix at
    Eigen::VectorXd snapshot_points(n_snapshots);  // a values in first row, b values in second row
    double *seq_value;
    for (int i = 0; i < n_snapshots; i++)
    {
        seq_value = halton(i, 1);
        snapshot_points(i) = seq_value[0] * (b_range.second - b_range.first) + b_range.first;
    }
    delete [] seq_value;

    // Compute FOM solutions and add them to the snapshot matrix
    for (int i = 0; i < n_snapshots; i++)
    {
        double b = snapshot_points(i);
        std::vector<double> solution = fom_solver.solve(b, t_eval);
        snapshot_matrix.conservativeResize(solution.size(), snapshot_matrix.cols() + 1);
        
        for (int j=0; j < n_snapshots; j++)
        {
            snapshot_matrix(j, i) = solution[j];
        }
    }
}

void Snapshots::compute_pod_basis()
{
    reference_state = snapshot_matrix.rowwise().mean();
    Eigen::MatrixXd centered_snapshots = snapshot_matrix.colwise() - reference_state;
    Eigen::BDCSVD<Eigen::MatrixXd, Eigen::DecompositionOptions::ComputeThinU> svd(centered_snapshots);
    basis = svd.matrixU();
}

void write_matrix(std::string &fileName, Eigen::MatrixXd  &matrix)
{
    const static Eigen::IOFormat csv_format(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
    std::ofstream file(fileName);
    if (file.is_open())
    {
        file << matrix.format(csv_format);
        file.close();
    }
}

Eigen::MatrixXd Snapshots::get_snapshot_matrix() {return snapshot_matrix;}

Eigen::MatrixXd Snapshots::get_basis() {return basis;}

Eigen::VectorXd Snapshots::get_reference_state() {return reference_state;}








