#include <iostream>
#include <cmath>
#include <vector>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

const double pi = 3.14159265358979323846;

double f(double x) {
    return std::sin(x);
}

int main() {
    int num_elements = 100;
    int num_nodes = num_elements + 1;
    double dx = (2 * pi) / num_elements;

    std::vector<std::vector<double>> K(num_nodes, std::vector<double>(num_nodes, 0));
    std::vector<double> F(num_nodes, 0);

    // Assemble the global stiffness matrix and load vector
    for (int i = 0; i < num_elements; ++i) {
        double x = i * dx;
        double x_next = (i + 1) * dx;

        K[i][i] += (1 + dx / 2) / dx;
        K[i][i + 1] += (-1 + dx / 2) / dx;
        K[i + 1][i] += (-1 - dx / 2) / dx;
        K[i + 1][i + 1] += (1 - dx / 2) / dx;

        F[i] += f(x) * dx / 2 + f(x_next) * dx / 2;
        F[i + 1] += f(x) * dx / 2 + f(x_next) * dx / 2;
    }

    // Modify the system to satisfy the u'(0) = 1 boundary condition using the method of Lagrange multipliers
    K[0][0] += 1;
    F[0] += dx;

    // Modify the system to satisfy the u(2*pi) = 0 boundary condition
    K[num_nodes - 1].assign(num_nodes, 0);
    K[num_nodes - 1][num_nodes - 1] = 1;
    F[num_nodes - 1] = 0;

    // Solve the system using Gaussian elimination
    for (int i = 0; i < num_nodes - 1; ++i) {
        for (int j = i + 1; j < num_nodes; ++j) {
            double factor = K[j][i] / K[i][i];
            for (int k = i; k < num_nodes; ++k) {
                K[j][k] -= factor * K[i][k];
            }
            F[j] -= factor * F[i];
        }
    }

    std::vector<double> u(num_nodes);
    u[num_nodes - 1] = F[num_nodes - 1] / K[num_nodes - 1][num_nodes - 1];
    for (int i = num_nodes - 2; i >= 0; --i) {
        double sum = 0;
        for (int j = i + 1; j < num_nodes; ++j) {
            sum += K[i][j] * u[j];
        }
        u[i] = (F[i] - sum) / K[i][i];
    }

    // Output the results and store x and u(x) values for plotting
    std::vector<double> x_values(num_nodes);
    std::vector<double> u_values(num_nodes);
    for (int i = 0; i < num_nodes; ++i) {
        double x = i * dx;
        std::cout << "Node " << i << ": x = " << x << ", u(x) = " << u[i] << std::endl;
        x_values[i] = x;
        u_values[i] = u[i];
    }

    // Plot the results
    plt::figure_size(800, 600);
    plt::plot(x_values, u_values, "b-o");
    plt::xlabel("x");
    plt::ylabel("u(x)");
    plt::title("FEM Solution for the Oscillatory Equation");
    plt::grid(true);
    plt::show();

    return 0;
}