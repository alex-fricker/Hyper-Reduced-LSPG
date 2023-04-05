#include "burgers_rewienski.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include "misc/eigen_utils.hpp"


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
    for (int i=0; i <= nx; i++)  // Adding extra cell for the ghost cell
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

Eigen::VectorXd BurgersRewienski::solve(float b) 
{
    std::cout << "\n-------------------\nStarting Solve\n-------------------\n" << std::endl;

    Eigen::VectorXd u0 = Eigen::VectorXd::Ones(nx + 1);  // Solution at pseudo time-step n
    Eigen::VectorXd u1 = Eigen::VectorXd::Ones(nx + 1);  // Solution at pseudo time-step n+1

    double time = 0;
    double dt;  //dt for current timestep

    residual = Eigen::VectorXd::Zero(nx);  // Residual in each cell for the current iteration
    jacobian = Eigen::MatrixXd::Ones(nx, nx);

    set_boundary_condition(u0, u1);
    set_initial_condition(u0);

    double tol = 1e-12;
    unsigned int max_itr = 5e4;
    unsigned int itr = 0;
    double residual_norm = 1;
    while (residual_norm > tol && itr < max_itr)
    {
        dt = set_timestep(u0);
        u1 = step_in_time(u0, b, dt);  // March solution forwards in time

        evaluate_residual(u1, residual, b);  // Compute Residual for current timestep
        residual_norm = residual.lpNorm<Eigen::Infinity>();  // L_inf norm of the residual
        normalized_residual.push_back(residual_norm);  // Normalize residual to track convergence

        if (itr % 1000 == 0)
        {
        std::cout << "\tIteration: " << itr
                  << "\n\tResidual: " << normalized_residual.back() << "\n"
                  << std::endl;
        }

        time += dt;
        itr++;
        u0 = u1;
    }
    std::cout << "-------------------\nDone Solving\n-------------------\n"
        << "Number of iterations: " << itr
        << "\nFinal non-linear residual: " << normalized_residual.back() 
        << std::endl;

    evaluate_jacobian(u1, jacobian, b);

    return u1(Eigen::seq(0, Eigen::placeholders::last - 1));
}

Eigen::VectorXd BurgersRewienski::step_in_time(const Eigen::VectorXd &u, const float &b, const double &dt)
{
    Eigen::VectorXd u1 = Eigen::VectorXd::Ones(u.size());
    double uhm;  // u_(n+0.5, j-0.5)
    double uhp;  // u_(n+0.5, j+0.5)

    for (int i = 1; i < nx; i++)
    {
        uhm = half_timestep_minus(u, i, dt, b);  // u_{i-0.5}^{n+0.5}
        uhp = half_timestep_plus(u, i, dt, b);  // u_{i+0.5}^{n+0.5}
        u1(i) = u(i) - dt / dx * (flux(uhp) - flux(uhm)) + dt * source_term(x[i], b);  // u_(n+1, j)
    }
    u1(nx) = u1(nx-1);  // Setting value of ghost cell
    u1(0) = u(0);
    
    return u1; 
}

double BurgersRewienski::half_timestep_plus(const Eigen::VectorXd &u, const unsigned int &i, const double &dt, const float &b) const
{
    return 0.5 * (u[i+1] + u[i]) - 0.5 * dt / dx * (flux(u[i+1]) - flux(u[i])) + 0.5 * dt * source_term(x[i] + 0.5 * dx, b);
}

double BurgersRewienski::half_timestep_minus(const Eigen::VectorXd &u, const unsigned int &i, const double &dt, const float &b) const
{
    return 0.5 * (u[i] + u[i-1]) - 0.5 * dt / dx * (flux(u[i]) - flux(u[i-1])) + 0.5 * dt * source_term(x[i] - 0.5 * dx, b);
}

void BurgersRewienski::evaluate_residual(
    const Eigen::VectorXd &u, 
    Eigen::VectorXd &residual,
    const float &b)
{
    double dt = set_timestep(u);
    for (int i = 1; i < nx - 1; i++)
    {
        double uhm = half_timestep_minus(u, i, dt, b);
        double uhp = half_timestep_plus(u, i, dt, b);
        double R = 1 / dx * (flux(uhp) - flux(uhm)) - source_term(x[i], b);
        residual(i) = R;
    }
}

void BurgersRewienski::evaluate_jacobian(
    const Eigen::VectorXd &u,
    Eigen::MatrixXd &jacobian,
    const double &b)
{
    double dt = set_timestep(u);
    for (int i = 1; i < nx - 1; i++)
    {
        jacobian(i, i) = (
            1 / dx * (0.5 * half_timestep_plus(u, i, dt, b) * (1 + dt / dx * u[i]) - 
                0.5 * half_timestep_minus(u, i, dt, b) * (1 - dt / dx * u[i])) 
        );
        if (i < nx - 1)
        {
            jacobian(i, i+1) = (
                -0.5 / dx * half_timestep_plus(u, i, dt, b) * (1 + dt / dx * u[i])
            );
        }
        if (i > 1)
        {
            jacobian(i, i-1) = (
                0.5 / dx * half_timestep_minus(u, i, dt, b) * (1 - dt / dx * u[i])
            );
        }
    }
}

void BurgersRewienski::set_boundary_condition(Eigen::VectorXd &u0, Eigen::VectorXd &u1)
{
	if (bc_spec.first == "Constant")  // Constant boundary condition
	{
		u0(0) = bc_spec.second[0];
        u1(0) = bc_spec.second[0];
        u0(nx) = bc_spec.second[0];
        u1(nx) = bc_spec.second[0];
	}
}

void BurgersRewienski::set_initial_condition(Eigen::VectorXd &u)
{
    if (ic_spec.first == "Constant")  // Constant initial condition
    {
        for (int i=0; i < nx; i++)
        {
            u(i) = ic_spec.second[0];
        }
    }
}

double BurgersRewienski::set_timestep(const Eigen::VectorXd &u)
{
    float dt = 1 * dx / u.lpNorm<Eigen::Infinity>();
    return dt;
}


double BurgersRewienski::flux(const double &u) const { return 0.5 * std::pow(u, 2); }

double BurgersRewienski::source_term(const double &x, const float &b) const {return 0.02 * std::exp(x * b); }

const Eigen::VectorXd BurgersRewienski::get_residual() const { return residual; }

const Eigen::MatrixXd BurgersRewienski::get_jacobian() const { return jacobian; }

const std::vector<double> BurgersRewienski::get_residual_history() const { return normalized_residual; }
























