#ifndef __BURGERS_REWIENSKI__
#define __BURGERS_REWIENSKI__

#include<vector>
#include<string>
#include "../Libraries/eigen/Eigen/Dense"

class BurgersRewienski
{
public:
	// Constructor
	BurgersRewienski(
        const int nx,
        const float x0,
        const float x1, 
		std::pair<std::string, std::vector<float>> bc_spec,
        std::pair<std::string, std::vector<float>> ic_spec);

	// Destructor
	~BurgersRewienski() {};

	// Function to evaluate solution
	std::vector<double> solve(float b, float t1);

	// Function to write solution to a text file
	void write_solution(const std::string &name, const std::vector<double> &u);

    const Eigen::VectorXd get_residual() const;

    const Eigen::MatrixXd get_jacobian() const;

    Eigen::VectorXd rom_residual(const Eigen::VectorXd &u, const float &b);

    Eigen::MatrixXd rom_jacobian(const Eigen::VectorXd &u, const float &b);

    // Constructs the residual vector for the current iteration
    void evaluate_residual(
        const  std::vector<std::vector<double>> &u, 
        Eigen::VectorXd &residual,  
        const double &dt,
        const float &b);

    // // Constructs the jacobian matrix dR/du for the current final solution
    // void evaluate_jacobian(
    //     const std::vector<std::vector<double>> &u,
    //     Eigen::MatrixXd &jacobian,
    //     const double &dt,
    //     const double &b);

    // Evaluate Burgers Flux
    double flux(const double &u);

    // Evaluate Burgers source term
    double source_term(const double &x, const float &b) const;

    // Number of elements
    int nx;

    // Boundary condition
    std::pair<std::string, std::vector<float>> bc_spec;

    // Initial condition
    std::pair<std::string, std::vector<float>> ic_spec;

    // Gridsize
    double dx;

    // Normalized residual at each timestep
    std::vector<double> normalized_residual;


private:
	// Function to set the boundary condition
	void set_boundary_condition(std::vector<std::vector<double>> &u);

	// Fucntion to set the initial condition
	void set_initial_condition(std::vector<std::vector<double>> &u);

    // Function to step forwards in time
    std::vector<double> step_in_time(const std::vector<double> &u, const float &b, const double &dt);

    // Determine timestep for the current iteration
    double set_timestep(const std::vector<double> &u);

    // Gridpoints
    std::vector<double> x;

    // Residual vector at the last iteration of the solution
    Eigen::VectorXd solution_residual;

    // Jacobian of the current solution
    Eigen::MatrixXd jacobian;

};
#endif 
