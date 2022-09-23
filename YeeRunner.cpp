#include <limits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/time.h>

// #include <Kokkos_Core.hpp>
#include <cmath>
#include "YeeGrid.h"

void checkSizes(int &N, int &M, int &S, int &nrepeat);

int main(int argc, char *argv[])
{
    std::cout << "foo";
    // Physical constants
    const double M_electron = 9.109e-31; // [kg]
    const double Q_electron = 1.602e-19; // [C] (intentionally positive)
    const double c = 2.99792458e8;
    const double mu_0 = 1.25663706e-6;
    const double eps_0 = 8.854187817e-12;
    const double k_B = 1.380649e-23; // Boltzmann constant

    // Nondimensionalized units
    const double M = M_electron; // [kg]
    const double Q = Q_electron; // [C]

    const double kappa = 5e1;                                                      // Nondimensional wave speed [m/s]
    const double lambda_D = 1e-2;                                                  // Debye length [m]
    const double n_bar = (M * eps_0 * pow(c, 2)) / pow((kappa * Q * lambda_D), 2); // Average macroscopic number density [m^-3]
    const double w_p = sqrt((n_bar * pow(Q, 2)) / (M * eps_0));                    // angular frequency
    const double T_bar = (n_bar * pow(Q * lambda_D, 2)) / (eps_0 * k_B);

    // More nondimensionalized units
    const double L = lambda_D; // In meters [m]
    const double T = 1 / w_p;  // In seconds/radians [s/r]
    const double V = L / T;    // In [m/s] (thermal velocity lam_D*w_p)

    // nondimensionalization parameters for involutions
    const double sig_1 = (M * eps_0) / (pow(Q, 2) * pow(T, 2) * n_bar);
    const double sig_2 = mu_0 * pow(Q, 2) * pow(L, 2) * n_bar / M;

    const double q_elec = -Q_electron / Q;
    const double q_ion = Q_electron / Q;
    const double m_elec = M_electron / M;
    const double m_ion = 1836 * m_elec;

    const char eleTag = '-';
    const char ionTag = '+';

    const double T_final = 10; // normalized wrt 1/w_p (plasma period)

    const int numGridIterations = 5;

    int grids[numGridIterations] = {8,16,32,64,128};
    double heatList[2][numGridIterations];
    double dxs[numGridIterations] = {0};

    // Physical grid parameters
    const double L_x = 50;
    const double L_y = 50;

    const double a_x = -25;
    const double b_x = a_x + L_x;

    const double a_y = -25;
    const double b_y = a_y + L_y;
    // End physical grid parameters

    for (int g = 0; g < numGridIterations; g++)
    {
        // Set up nondimensional grid
        const int Nx = grids[g];
        const int Ny = grids[g];
        const double dx = (b_x-a_x)/(Nx-1);
        const double dy = (b_y-a_y)/(Ny-1);

        double x[Nx];
        double y[Ny];
        for (int i = 0; i < Nx; i++) {
            x[i] = a_x + i*dx;
            y[i] = a_y + i*dy;
        }
        double Hz[Nx][Ny];
        double Ex[Nx][Ny];
        double Ey[Nx][Ny];
        // End set up

        dxs[g] = dx;
        double dt = 1 / (2 * kappa * sqrt(1 / pow(dx, 2) + 1 / pow(dy, 2)));

        int N_steps = floor(T_final / dt);

        const double Jx_stagger_x = dx / 2;
        const double Jx_stagger_y = 0;
        const double Jy_stagger_x = 0;
        const double Jy_stagger_y = dy / 2;

        const int Jx_padding_x = 1;
        const int Jx_padding_y = 0;
        const int Jy_padding_x = 0;
        const int Jy_padding_y = 1;

        const double Hz_stagger_x = dx / 2;
        const double Hz_stagger_y = dy / 2;

        const double offset_x = -a_x;
        const double offset_y = -a_y;

        YeeGrid yeeGrid(Nx, Ny, dx, dy, dt, kappa);

        for (int n = 0; n < N_steps; n++) {
            if (n % 50 == 0) {
                std::cout << std::to_string(n) << "\n";
                yeeGrid.print();
            }
            yeeGrid.step();
        }
        return 0;
    }
        
    return 0;
}
