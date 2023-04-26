#include <iostream>
#include <vector>
#include <cmath>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
using namespace std;

// Analytic solution function
double analytic_solution(double x) {
    return (3/(2*exp(2*M_PI)) + 0.5 - 3*exp(-x)/2 + (-cos(x)-sin(x))/2);
}

// Left and right bounds
const double left_bound = 0;
const double right_bound = 2 * M_PI;

// Number of elements and length of each element
int num_elements = 40;
double element_length = (right_bound - left_bound) / num_elements;

// Number of nodes (equal to number of elements plus one)
int num_nodes = num_elements + 1;

// Left-hand side (LHS) matrix
vector<vector<double>> LHS(num_nodes, vector<double>(num_nodes));

// Right-hand side (RHS) vector
vector<double> RHS(num_nodes, 0.0);

// Assemble the LHS matrix
void assemble_LHS() {
    for (int i = 0; i < num_elements; i++) {
        double inv_delta_x = 1.0 / element_length;
        vector<vector<double>> element_matrix = {{-inv_delta_x - 0.5, inv_delta_x + 0.5}, {inv_delta_x - 0.5, -inv_delta_x + 0.5}};
        LHS[i][i] += element_matrix[0][0];
        LHS[i][i + 1] += element_matrix[0][1];
        LHS[i + 1][i] += element_matrix[1][0];
        LHS[i + 1][i + 1] += element_matrix[1][1];
    }
}

// Left integral 1
double left_integral_1(double x, double xu, double inv_delta_x) {
    return inv_delta_x * (xu * (-cos(x)) + x * cos(x) - sin(x));
}

// Left integral 2
double left_integral_2(double x, double xl, double inv_delta_x) {
    return inv_delta_x * (-x * cos(x) + sin(x) + xl * cos(x));
}

// Assemble the RHS vector
void assemble_RHS() {
    for (int i = 0; i < num_elements; i++) {
        double xu = (i + 1) * element_length;
        double xl = i * element_length;
        double delta_x = xu - xl;
        double inv_delta_x = 1.0 / element_length;
        vector<double> element_rhs = {left_integral_1(xu, xu, inv_delta_x) - left_integral_1(xl, xu, inv_delta_x), left_integral_2(xu, xl, inv_delta_x) - left_integral_2(xl, xl, inv_delta_x)};
        RHS[i] += element_rhs[0];
        RHS[i + 1] += element_rhs[1];
    }

    // Apply the right boundary condition
    double right_boundary_condition = 0.0;
    RHS[num_elements] = right_boundary_condition;
    LHS[num_elements][num_elements - 1] = 0.0;
    LHS[num_elements - 1][num_elements] = 0.0;
    LHS[num_elements][num_elements] = 1.0;
    RHS[0] = RHS[0] + 1.0;
}

// Solve the system of equations and compute the solution at each node
void solve(std::vector<double>& x_values, std::vector<double>& u_values) {
    // Create a vector to hold the solution
    vector<double> solution(num_nodes);
    for (int i = 0; i < num_nodes; i++) {
        solution[i] = 0.0;
    }

    // Eliminate variables using Gaussian elimination with partial pivoting
    for (int k = 0; k < num_nodes - 1; k++) {
        for (int i = k + 1; i < num_nodes; i++) {
            double factor = LHS[i][k] / LHS[k][k];
            for (int j = k + 1; j < num_nodes; j++) {
                LHS[i][j] -= factor * LHS[k][j];
            }
            RHS[i] -= factor * RHS[k];
        }
    }

    // Back-substitute to find the solution
    for (int i = num_nodes - 1; i >= 0; i--) {
        solution[i] = RHS[i] / LHS[i][i];
        for (int j = i - 1; j >= 0; j--) {
            RHS[j] -= LHS[j][i] * solution[i];
        }
    }

    // Collect the data points for plotting
    x_values.clear();
    u_values.clear();
    x_values.reserve(num_nodes);
    u_values.reserve(num_nodes);

    for (int i = 0; i < num_nodes; i++) {
        double x = left_bound + i * element_length;
        x_values.push_back(x);
        u_values.push_back(solution[i]);
        cout << "u[" << x << "] = " << solution[i] << endl;
    }
}

int main() {
    // Assemble the LHS matrix and RHS vector for 40 elements
    assemble_LHS();
    assemble_RHS();

    // Solve the system of equations and compute the solution at each node for 40 elements
    cout << num_elements << " Elements Solution: " << endl;
    std::vector<double> x_values1, u_values1;
    solve(x_values1, u_values1);

    // Assemble the LHS matrix and RHS vector for 20 elements
    num_elements = 20;
    element_length = (right_bound - left_bound) / num_elements;
    num_nodes = num_elements + 1;
    LHS = vector<vector<double>>(num_nodes, vector<double>(num_nodes));
    RHS = vector<double>(num_nodes, 0.0);
    assemble_LHS();
    assemble_RHS();

    // Solve the system of equations and compute the solution at each node for 20 elements
    cout << num_elements << " Elements Solution: " << endl;
    std::vector<double> x_values2, u_values2;
    solve(x_values2, u_values2);

    // Plot the results
    plt::figure_size(1200, 800);
    plt::named_plot("40 Equal Steps", x_values1, u_values1, "b.-");
    plt::named_plot("20 Elements", x_values2, u_values2, "r.-");
    plt::xlabel("x");
    plt::ylabel("u(x)");
    plt::title("Finite Element Method Solution");
    plt::legend();
    plt::show();

    // Plot the analytical solution
    int num_points = 100;
    double dx = (2 * M_PI) / (num_points - 1);
    std::vector<double> x_values(num_points), y_values(num_points);
    for (int i = 0; i < num_points; ++i) {
        double x = i * dx;
        x_values[i] = x;
        y_values[i] = analytic_solution(x);
    }

    plt::figure_size(1200, 800);
    plt::xlabel("x");
    plt::ylabel("u(x)");
    plt::title("Analytic Solution");
    plt::named_plot("Analytic Solution", x_values, y_values);
    plt::show();

    return 0;
}