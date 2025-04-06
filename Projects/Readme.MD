# HJI
This directory contains an implementation of the Bellman equation from the paper:
"A Minimum Discounted Reward Hamilton–Jacobi Formulation for Computing Reachable Sets",
using the value iteration method.

It includes two Julia files related to:

Hamilton-Jacobi-Bellman (HJB): A special case of HJI without disturbance (i.e., no second player involved).

Hamilton-Jacobi-Isaacs (HJI): A more general formulation that solves the variational inequality for the HJI partial differential equation, accounting for disturbances or adversarial inputs (second player).

# Average Reward (AR) MDP
This directory contains two subdirectories: LP_primal and LP_dual.
AR MDP problem is formulated as a linear programming (LP) problem. Using JuMP and Gurobi, both the primal and dual formulations of the optimization problem are implemented and solved.

# General Notes
The system under consideration is a simple double integrator with disturbance affecting the second state. The disturbance is modeled as a normally distributed random variable.

The sampling resolution is set to 100. In the Average Reward MDP cases, when a state lies between grid points, the nearest neighbor is selected using the K-Nearest Neighbors (KNN) method.

In contrast, for the HJI implementation, interpolation is used to estimate values between grid points. This is done via a linear combination of the four nearest grid points, weighted by their respective probabilities.


