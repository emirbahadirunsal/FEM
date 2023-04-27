# FEM

# Homework 1

This repository contains the code for solving the following oscillatory ODE with boundary conditions:

    ```
    d^2u/dx^2 + du/dx = sin(x),
    u(2*pi) = 0,    ```
    u'(0) = 1
    ```

The analytical solution of the ODE is given by:

    ```
    u(x) = 3/(2*e^(2*pi)) + 1/2 - 3/2*e^(-x) - (cos(x) + sin(x))/2
    ```

## Files

### `homework1.cpp`

This C++ script uses the finite difference method to solve the differential equation and plot output the solution to the console.

### `homework1_plot.cpp`

This will generate the solution plots.

### `fem_solution.eps`

FEM result of our code.

### `analytic_solution.eps`

Analytical result for comparison.

## Usage

Run the following command to compile `homework1_plot.cpp`:
   ```
   g++ homework1_plot.cpp -std=c++11 -I/usr/local/include/python2.7 -I/usr/local/lib/python2.7/site-packages/numpy/core/include -lpython2.7
   ```

Note that the programs require the following libraries to be installed:
- `matplotlibcpp`: A C++ library for plotting using Matplotlib in Python.
- `numpy`: A Python library for numerical computing.
