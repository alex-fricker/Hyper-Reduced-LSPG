#include "burgers_rewienski.hpp"
#include <iostream>
#include <fstream>
#include <cmath>

BurgersRewienski::BurgersRewienski(
		const int nx, 
        const float x0, 
        const float x1,
        std::pair<std::string, std::vector<float>> bc_spec,
        std::pair<std::string, std::vector<float>> ic_spec)
	: nx(nx)
    , bc_spec(bc_spec)
    , ic_spec(ic_spec)
{
	dx = (x1 - x0) / (nx - 1);
    for (int i=0; i <= nx; i++)
    {
		x.push_back(x0 + i * dx);
	}

    std::cout <<
        "Initializing Burgers Rewienski class with:\n\tSpacial domain " << x0 << " < x < " << x1
        << "\n\tnx = " << nx
        << "\n\tdx = " << dx
        << "\n\tBC type " << bc_spec.first
        << "\n\tIC type " << ic_spec.first
        << "\n\n";

}

std::vector<double> BurgersRewienski::solve(float b, float t1, bool return_residual) 
{
    std::cout << "\n-------------------\nStarting Solve\n-------------------\nEvaluating at t=" << t1 << std::endl;

    std::vector<std::vector<double>> u(2, std::vector<double>(x.size(), 1));  // Computational grid
    double time = 0;
    double dt = 0.5 * t1;

    set_boundary_condition(u);
    set_initial_condition(u);

    double time = 0;
    while (t1 - time > 0.0001)
    {
        double dt = set_timestep(u[0]);

        u[1] = step_in_time(u[0], b, dt);
        u[0] = u[1];

        std::cout << "\ttime: " << time 
            << "\n\tdt: " << dt << "\n" <<std::endl;

        time += dt;
    }

    // Compute normalized residual
    if (return_residual)
    {
        double R = 0;
        for (int i = 0; i < nx; i++)
        {
            R += std::pow(2, cell_residual(u[1][i], u[0][i], u[1][i+1], u[1][i-1], x[i], b, dt)) * dx;
        }
        normalized_residual = std::pow(0.5, R);
    }

    std::cout << "Done Solving\n" << std::endl;
    return u[1];
}

void BurgersRewienski::set_boundary_condition(std::vector<std::vector<double>> &u)
{
	if (bc_spec.first == "Constant")  // Constant boundary condition
	{
		u[0][0] = bc_spec.second[0];
        u[1][0] = bc_spec.second[0];
        u[0][u[0].size()-1] = bc_spec.second[0];
        u[1][u[0].size()-1] = bc_spec.second[0];
	}
}

void BurgersRewienski::set_initial_condition(std::vector<std::vector<double>> &u)
{
    if (ic_spec.first == "Constant")  // Constant initial condition
    {
        for (int i=0; i < nx; i++)
        {
            u[0][i] = ic_spec.second[0];
        }
    }
}

std::vector<double> BurgersRewienski::step_in_time(const std::vector<double> &u, const float &b, const double &dt)
{
    std::vector<double> u1(u.size(), 1);
    double uhm;  // u_(n+0.5, j-0.5)
    double uhp;  // u_(n+0.5, j+0.5)

    for (int i = 1; i < nx; i++)
    {
        uhm = 0.5 * (u[i] + u[i-1]) - 0.5 * dt / dx * (flux(u[i]) - flux(u[i-1])) + 0.5 * dt * source_term(x[i] - 0.5 * dx, b);
        uhp = 0.5 * (u[i+1] + u[i]) - 0.5 * dt / dx * (flux(u[i+1]) - flux(u[i])) + 0.5 * dt * source_term(x[i] + 0.5 * dx, b);
        u1[i] = u[i] - dt / dx * (flux(uhp) - flux(uhm)) + dt * source_term(x[i], b);  // u_(n+1, j)
    }
    u1[u1.size() - 1] = u1[u1.size() - 2];  // Setting value of ghost cell
    u1[0] = u[0];
    
    return u1; 
}

double BurgersRewienski::set_timestep(const std::vector<double> &u)
{
    float u_max = 0;
    for (int i=0; i < nx; i++)
    {
        u_max = std::abs(u[i]) > u_max ? std::abs(u[i]) : u_max;
    }
    float dt = 0.999999 * dx / u_max;
    return dt;
}

void BurgersRewienski::write_solution(const std::string &name, const std::vector<double> &u)
{
    std::ofstream file;
    std::string fname = name + ".txt";
    file.open(fname );
    file << "Value, x position\n";
    for (int i=0; i < nx; i++)
    {
        file << u[i] << ',' << x[i] << "\n";
    }
    file.close();
}

double BurgersRewienski::flux(const double &u) { return 0.5 * std::pow(u, 2); }

double BurgersRewienski::source_term(const double &x, const float &b) {return 0.02 * std::exp(x * b); }

double BurgersRewienski::cell_residual(
    const double &u11, 
    const double &u10, 
    const double &u21,
    const double &u01,
    const double &x,
    const float &b, 
    const double dt)
{
    double R = (u11 - u10) / dt + u11 * (u21 - u01) / (2 * dx) - source_term(x, b);
    return R;
}

double BurgersRewienski::get_residual()
{
    return normalized_residual;
}
























