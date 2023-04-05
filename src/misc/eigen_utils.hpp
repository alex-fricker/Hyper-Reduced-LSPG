#ifndef __EIGEN_UTILS__
#define __EIGEN_UTILS__

#include <vector>
#include "../../Libraries/eigen/Eigen/Dense"

void write_matrix(const std::string &fileName, const Eigen::MatrixXd &matrix);
void write_vector(const std::string &fileName, const Eigen::VectorXd &matrix);
std::vector<double> eigen2std_vector(const Eigen::VectorXd &v);
Eigen::VectorXd std2eigen_vector(const std::vector<double> &v);

#endif
