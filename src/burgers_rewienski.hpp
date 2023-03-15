#ifndef __BURGERS_REWIENSKI__
#define __BURGERS_REWIENSKI__

#include<vector>
#include<string>

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

    const std::vector<double> get_residual() const;

    // Gridsize
    double dx;

    // Number of elements
    int nx;

private:
	// Function to set the boundary condition
	void set_boundary_condition(std::vector<std::vector<double>> &u);

	// Fucntion to set the initial condition
	void set_initial_condition(std::vector<std::vector<double>> &u);

    // Function to step forwards in time
    std::vector<double> step_in_time(const std::vector<double> &u, const float &b, const double &dt);

    // Evaluate Burgers Flux
    double flux(const double &u);

    // Evaluate Burgers source term
    double source_term(const double &x, const float &b) const;

    // Determine timestep for the current iteration
    double set_timestep(const std::vector<double> &u);

    // Evaluate the full order residual
    double cell_residual(
        const double &u11, 
        const double &u10, 
        const double &u21,
        const double &u01,
        const double &x,
        const float &b, 
        const double dt) const;

    // Gridpoints
    std::vector<double> x;

    // Boundary condition
    std::pair<std::string, std::vector<float>> bc_spec;

    // Initial condition
    std::pair<std::string, std::vector<float>> ic_spec;

    // Normalized residual at each timestep
    std::vector<double> normalized_residual;
};
#endif 
