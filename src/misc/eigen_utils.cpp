#include <iostream>
#include <fstream>
#include <vector>
#include "eigen_utils.hpp"

void write_matrix(const std::string &fileName, const Eigen::MatrixXd &matrix)
{
    const static Eigen::IOFormat csv_format(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
    std::ofstream file(fileName);
    if (file.is_open())
    {
        file << matrix.format(csv_format);
        file.close();
    }
}

void write_vector(const std::string &fileName, const Eigen::VectorXd &matrix)
{
    const static Eigen::IOFormat csv_format(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
    std::ofstream file(fileName);
    if (file.is_open())
    {
        file << matrix.format(csv_format);
        file.close();
    }
}

std::vector<double> eigen2std_vector(const Eigen::VectorXd &v)
{
    std::vector<double> vec;
    for (int i = 0; i < v.size(); i++)
    {
        vec.push_back(v(i));
    }
    return vec;
}

Eigen::VectorXd std2eigen_vector(const std::vector<double> &v)
{
    Eigen::VectorXd vec(v.size());
    int size = static_cast<int>(v.size());
    for (int i = 0; i < size; i++)
    {
        vec(i) = v[i];
    }
    return vec;
}
