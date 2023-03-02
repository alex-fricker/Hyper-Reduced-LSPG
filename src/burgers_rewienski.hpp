#ifndef __BURGERS_REWIENSKI__
#define __BURGERS_REWIENSKI__

#include<vector>
#include<string>

class BurgersRewienski
{
public:
	// Constructor
	BurgersRewienski(
        const int nx, const float x0, const float x1, const float CFL, 
		std::pair<std::string, std::vector<float>> bc_spec,
        std::pair<std::string, std::vector<float>> ic_spec);

	// Destructor
	~BurgersRewienski() {};

	// Function to evaluate solution
	std::vector<double> solve(float b, float t1);

	// Function to write solution to a text file
	void write_solution(const std::string &name, const std::vector<double> &u);

private:
	// Function to set the boundary condition
	void set_boundary_condition(std::vector<std::vector<double>> &u);

	// Fucntion to set the initial condition
	void set_initial_condition(std::vector<std::vector<double>> &u);

    // Function to set the timestep for the current iteration
    double set_timestep(const std::vector<double> &u) const;

    // Function to step forwards in time
    std::vector<double> step_in_time(const std::vector<double> &u, const float &b, const double &dt);

    // Gridpoints
    std::vector<double> x;

    // Gridsize
    double dx;

    // Number of elements
    int nx;

    // CFL
    float CFL;

    // Boundary condition
    std::pair<std::string, std::vector<float>> bc_spec;

    // Initial condition
    std::pair<std::string, std::vector<float>> ic_spec;

};
#endif 
