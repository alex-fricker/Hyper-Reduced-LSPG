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
	Eigen::VectorXd solve(float b);

    const Eigen::VectorXd get_residual() const;

    const Eigen::MatrixXd get_jacobian() const;

    const std::vector<double> get_residual_history() const;

    // Constructs the steady state residual vector for the given iteration
    void evaluate_residual(
        const  Eigen::VectorXd &u, 
        Eigen::VectorXd &residual,
        const float &b);

    // Constructs the jacobian matrix dR/du for the current final solution
    void evaluate_jacobian(
        const Eigen::VectorXd &u,
        Eigen::MatrixXd &jacobian,
        const double &b);

    // Evaluate Burgers Flux
    double flux(const double &u) const;

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


private:
	// Function to set the boundary condition
	void set_boundary_condition(Eigen::VectorXd &u0, Eigen::VectorXd &u1);

	// Fucntion to set the initial condition
	void set_initial_condition(Eigen::VectorXd &u);

    // Function to step forwards in time
    Eigen::VectorXd step_in_time(const Eigen::VectorXd &u, const float &b, const double &dt);

    // Determine timestep for the current iteration
    double set_timestep(const Eigen::VectorXd &u);

    // Step 1 for the two step Lax-Wendroff, evaluates u_{i+0.5}^{n+0.5}
    double half_timestep_plus(const Eigen::VectorXd &u, const unsigned int &i, const double &dt, const float &b) const;

    // Step 1 for the two step Lax-Wendroff, evaluates u_{i-0.5}^{n+0.5}
    double half_timestep_minus(const Eigen::VectorXd &u, const unsigned int &i, const double &dt, const float &b) const;

    // Gridpoints
    std::vector<double> x;

    // Residual vector at the last iteration of the solution
    Eigen::VectorXd residual;

    // Normalized residual at each pseudo-timestep
    std::vector<double> normalized_residual;

    // Jacobian of the current solution
    Eigen::MatrixXd jacobian;

};
#endif 
