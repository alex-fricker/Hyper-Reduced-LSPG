#include "burgers_rewienski.hpp"
#include <iostream>
#include <fstream>
#include <cmath>

BurgersRewienski::BurgersRewienski(
		const int nx, const float x0, const float x1, const float CFL, 
        std::pair<std::string, std::vector<float>> bc_spec,
        std::pair<std::string, std::vector<float>> ic_spec)
	: nx(nx)
    , CFL(CFL)
    , bc_spec(bc_spec)
    , ic_spec(ic_spec)
{
	dx = (x1 - x0) / (nx - 1);
	for (int i=0; i < nx; i++)
    {
		x.push_back(x0 + i * dx);
	}

    std::cout <<
        "Initializing Burgers Rewienski class with:\n\tSpacial domain " << x0 << " < x < " << x1
        << "\n\tnx = " << nx
        << "\n\tdx = " << dx
        << "\n\tCFL = " << CFL
        << "\n\tBC type " << bc_spec.first
        << "\n\tIC type " << ic_spec.first
        << "\n\tCFD: " << CFL << std::endl;

}

std::vector<float> BurgersRewienski::solve(float b, float t1) 
{
    std::cout << "\n-------------------\nStarting Solve\n-------------------\nEvaluating at t=" << t1 << std::endl;

    std::vector<std::vector<float>> u(2, std::vector<float>(x.size(), 1));  // Computational grid
    std::vector<float> t;  // Keep track of timesteps
    float time = 0;

    set_boundary_condition(u);
    set_initial_condition(u);

    while (time < t1)
    {
//        float dt = set_timestep(u[0]);
        float dt = 0.01;
        t.push_back(dt);
        if (time + dt > t1) dt = t1 - time;

        u[1] = step_in_time(u[0], b, dt);
        u[0] = u[1];

        std::cout << "\ttime: " << time 
            << "\n\tdt: " << dt << "\n" <<std::endl;

        time += dt;
    }
    std::cout << "Done Solving\n" << std::endl;
    return u[1];
}

void BurgersRewienski::set_boundary_condition(std::vector<std::vector<float>> &u)
{
	if (bc_spec.first == "Constant")  // Constant boundary condition
	{
		u[0][0] = bc_spec.second[0];
        u[1][0] = bc_spec.second[0];
        u[0][u[0].size()-1] = bc_spec.second[0];
        u[1][u[0].size()-1] = bc_spec.second[0];
	}
}

void BurgersRewienski::set_initial_condition(std::vector<std::vector<float>> &u)
{
    if (ic_spec.first == "Constant")  // Constant initial condition
    {
        for (int i=0; i < nx; i++)
        {
            u[0][i] = ic_spec.second[0];
        }
    }
}


float BurgersRewienski::set_timestep(const std::vector<float> &u) const
{
    float u_max = 0;
    for (int i=0; i < nx; i++)
    {
        u_max = std::abs(u[i]) > u_max ? std::abs(u[i]) : u_max;
    }
    float dt = CFL * dx / u_max;
    return dt;
}

std::vector<float> BurgersRewienski::step_in_time(const std::vector<float> &u, const float &b, const float &dt)
{
    std::vector<float> u_next = u;
    for (int i=1; i < nx; i++)
    {
        u_next[i] = (
            u[i] + 0.02 * exp(x[i] * b) - dt / (4 * dx) * 
            (
                pow(u[i+1], 2) - pow(u[i], 2)
            ) + pow(dt, 2) / (8 * pow(dx, 2)) * 
            (
                (u[i+1] + u[i]) * (pow(u[i+1], 2) - pow(u[i], 2)) - 
                (u[i] + u[i-1]) * (pow(u[i], 2) - pow(u[i-1], 2))
            )
        );
    }
    return u_next; 
}

void BurgersRewienski::write_solution(const std::string &name, const std::vector<float> &u)
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




























