#include "../../Libraries/eigen/Eigen/Dense"
#include <iostream>
#include <fstream>

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
