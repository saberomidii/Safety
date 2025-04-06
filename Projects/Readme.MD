# HJI
This directory contains an implementation of the Bellman equation from the paper:
"A Minimum Discounted Reward Hamilton–Jacobi Formulation for Computing Reachable Sets",
using the value iteration method.

It includes two Julia files related to:

Hamilton-Jacobi-Bellman (HJB): A special case of HJI without disturbance (i.e., no second player involved).

Hamilton-Jacobi-Isaacs (HJI): A more general formulation that solves the variational inequality for the HJI partial differential equation, accounting for disturbances or adversarial inputs (second player).

