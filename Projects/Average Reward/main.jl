##########################################################
# Average Reward and Linear Programming for Double-Integrator Systems
##########################################################
using Distributions
using StatsFuns
using Statistics
using Random
using NearestNeighbors
using JuMP
import Gurobi
using FileIO
using Base: mkpath

##########################################################
# 1) Setup: parameters, discretization, etc.
##########################################################

# Set Gurobi WLS license credentials
ENV["GRB_WLSACCESSID"] = "52eb20bf-115c-42d3-931f-47561460611c"
ENV["GRB_WLSSECRET"] = "587b8f93-6d53-43c9-af49-6b96ac589004"
ENV["GRB_LICENSEID"] = "2611020"

# Optional: set a random seed if desired
# Random.seed!(2)

# System parameters
dt = 0.1  # time step
g = 9.81  # gravity
m = 1.0   # mass
l = 1.0   # length

# State space discretization
nx = 50   # number of grid points in x-dimension
ny = 50   # number of grid points in y-dimension
xmin, xmax = -5.0, 5.0
ymin, ymax = -5.0, 5.0
x_grid = range(xmin, xmax, length=nx)
y_grid = range(ymin, ymax, length=ny)

# Action space discretization
na = 5    # number of discrete actions
amin, amax = -10.0, 10.0
actions = range(amin, amax, length=na)

# Noise parameters
noise_std = 0.1  # standard deviation of noise

# Create state-action pairs
states = [(x, y) for x in x_grid for y in y_grid]
sa_pairs = [(s, a) for s in states for a in actions]
num_states = length(states)
num_actions = length(actions)

##########################################################
# 2) Dynamics and reward functions
##########################################################

# Dynamics function - double integrator
function dynamics(s, a, w)
    x, y = s
    
    # Apply control and noise
    x_next = x + dt * y
    y_next = y + dt * (a + w) / m
    
    # Apply bounds
    x_next = clamp(x_next, xmin, xmax)
    y_next = clamp(y_next, ymin, ymax)
    
    return (x_next, y_next)
end

# Reward function
function reward(s, a)
    x, y = s
    # Negative quadratic cost on state and action
    return -(x^2 + 0.1*y^2 + 0.01*a^2)
end

##########################################################
# 3) Transition probability calculation
##########################################################

# Normal distribution for noise
noise_dist = Normal(0, noise_std)

# Function to find nearest state in the grid
function find_nearest_state(s)
    # Create KD-tree from states array for efficient nearest neighbor search
    kdtree = KDTree(reshape(reinterpret(Float64, states), (2, num_states)))
    
    # Find the index of the nearest state
    idx, dist = nn(kdtree, collect(Float64, s))
    
    return states[idx[1]]
end

# Calculate transition probabilities
function calculate_transition_probs()
    # Number of samples for Monte Carlo
    num_samples = 100
    
    # Initialize transition probability matrix
    P = zeros(num_states, num_states, num_actions)
    
    # Create mapping from state to index
    state_to_idx = Dict(state => i for (i, state) in enumerate(states))
    
    # For each state-action pair
    for (s_idx, s) in enumerate(states)
        for (a_idx, a) in enumerate(actions)
            # Monte Carlo sampling for stochastic transitions
            next_states = []
            for _ in 1:num_samples
                # Sample noise
                w = rand(noise_dist)
                
                # Apply dynamics
                s_next = dynamics(s, a, w)
                
                # Find nearest grid point
                s_next_grid = find_nearest_state(s_next)
                
                push!(next_states, s_next_grid)
            end
            
            # Count occurrences of each next state
            for s_next in next_states
                s_next_idx = state_to_idx[s_next]
                P[s_idx, s_next_idx, a_idx] += 1/num_samples
            end
        end
    end
    
    return P
end

# Calculate reward for each state-action pair
function calculate_rewards()
    R = zeros(num_states, num_actions)
    
    for (s_idx, s) in enumerate(states)
        for (a_idx, a) in enumerate(actions)
            R[s_idx, a_idx] = reward(s, a)
        end
    end
    
    return R
end

##########################################################
# 4) Linear Programming for Average Reward
##########################################################

function solve_average_reward()
    println("Calculating transition probabilities...")
    P = calculate_transition_probs()
    
    println("Calculating rewards...")
    R = calculate_rewards()
    
    println("Setting up linear program...")
    
    # Create JuMP model with Gurobi solver
    model = Model(Gurobi.Optimizer)
    
    # Variables: state-action occupancy measures
    @variable(model, x[1:num_states, 1:num_actions] >= 0)
    
    # Variable: average reward
    @variable(model, rho)
    
    # Objective: maximize average reward
    @objective(model, Max, rho)
    
    # Constraint: total probability = 1
    @constraint(model, sum(x) == 1)
    
    # Bellman flow constraints
    for j in 1:num_states
        @constraint(model, 
            sum(x[j, :]) == 
            sum(x[i, a] * P[i, j, a] for i in 1:num_states, a in 1:num_actions)
        )
    end
    
    # Average reward constraint
    @constraint(model, 
        rho == sum(x[s, a] * R[s, a] for s in 1:num_states, a in 1:num_actions)
    )
    
    # Solve the model
    println("Solving linear program...")
    optimize!(model)
    
    if termination_status(model) == MOI.OPTIMAL
        println("Optimal solution found!")
        
        # Extract optimal policy and value function
        x_opt = value.(x)
        rho_opt = value(rho)
        
        # Compute policy from occupancy measures
        policy = zeros(Int, num_states)
        for s in 1:num_states
            # For each state, find the action with highest occupancy measure
            _, policy[s] = findmax(x_opt[s, :])
        end
        
        return rho_opt, policy, x_opt
    else
        println("No optimal solution found. Status: ", termination_status(model))
        return nothing, nothing, nothing
    end
end

##########################################################
# 5) Main execution
##########################################################

function main()
    println("Starting Average Reward optimization for Double-Integrator System")
    
    # Create results directory
    results_dir = joinpath(@__DIR__, "results")
    mkpath(results_dir)
    
    # Solve the average reward problem
    rho, policy, occupancy = solve_average_reward()
    
    if !isnothing(rho)
        println("Optimal average reward: ", rho)
        
        # Save results
        open(joinpath(results_dir, "results.txt"), "w") do f
            println(f, "Optimal average reward: ", rho)
            println(f, "Policy:")
            for (s_idx, s) in enumerate(states)
                a_idx = policy[s_idx]
                println(f, "State: ", s, " → Action: ", actions[a_idx])
            end
        end
        
        println("Results saved to: ", joinpath(results_dir, "results.txt"))
    end
end

# Run the main function
main()
