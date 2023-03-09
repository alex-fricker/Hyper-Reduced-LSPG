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

std::vector<double> BurgersRewienski::solve(float b, float t1) 
{
    std::cout << "\n-------------------\nStarting Solve\n-------------------\nEvaluating at t=" << t1 << std::endl;

    std::vector<std::vector<double>> u(2, std::vector<double>(x.size(), 1));  // Computational grid
    // double time = 0;

    set_boundary_condition(u);
    set_initial_condition(u);

    // while (time < t1)
    // {
    //     double dt = dx;
    //     if (time + dt > t1) { dt = t1 - time; }

    //     u[1] = step_in_time(u[0], b, dt);
    //     u[0] = u[1];

    //     std::cout << "\ttime: " << time 
    //         << "\n\tdt: " << dt << "\n" <<std::endl;

    //     time += dt;
    // }
    double dt = t1;
    u[1] = step_in_time(u[0], b, dt);
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
    std::vector<double> u1(x.size(), 1);
    for (int i = 1; i < nx; i++)
    {
        u1[i] = LW_step2(u, x[i], i, dt, b);
    }
    u1[u1.size() - 1] = u[u.size() - 2];  // Setting value of ghost cell
    return u1; 
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

double BurgersRewienski::LW_step1(
    const std::vector<double> &u,
    const double &x,
    const int i,
    const double dt,
    const float b)
{
    double u1 = u[i] - dt / dx * (flux(u[i+1]) - flux(u[i])) + dt * source_term(x, b);
    return u1;
}


double BurgersRewienski::LW_step2(
    const std::vector<double> &u,
    const double &x,
    const int i,
    const double dt,
    const float b)
{
    double u1 = (
        u[i] - dt / dx * (flux(LW_step1(u, x, i, dt, b)) - flux(LW_step1(u, x, i, dt, b))) +
        dt * source_term(x, b)
    );
    return u1;
}

























