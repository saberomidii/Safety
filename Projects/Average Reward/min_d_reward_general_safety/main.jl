##########################################################
# Average Reward and Linear Programming for Double-Integrator Systems 
# Primal Formulation for Computing Reachable Sets
##########################################################
using Dates
using DelimitedFiles  # For saving data to CSV files
using Distributions
using StatsFuns
using Statistics
using Random
using NearestNeighbors
using JuMP
import Gurobi
using FileIO   # for saving results if needed
using Base: mkpath

##########################################################
# 1) Setup: parameters, discretization, etc.
##########################################################
# Set Gurobi WLS license credentials
println("Setting Gurobi license credentials...")
ENV["GRB_WLSACCESSID"] = "52eb20bf-115c-42d3-931f-47561460611c"
ENV["GRB_WLSSECRET"] = "587b8f93-6d53-43c9-af49-6b96ac589004"
ENV["GRB_LICENSEID"] = "2611020"
println("License credentials set")

# Optional: set a random seed if desired
# Random.seed!(50)

# Single noise standard deviation
const sigma = 1.0

const threshold_for_transit = 0.001

const x_1_min = -1.0
const x_1_max = 5.0
const x_2_min = -5.0
const x_2_max = 5.0

const u_min = -2.0
const u_max = 2.0

const d_min=-2.0
const d_max= 2.0 


const num_points_action = 11
const num_points_state = 161

# Number of random samples used when constructing transitions
const nsamples = 100

# Number of states in each dimension
x1 = collect(LinRange(x_1_min, x_1_max, num_points_state))  # Possible x values
x2 = collect(LinRange(x_2_min, x_2_max, num_points_state))  # Possible v values

# Cartesian product of (x, v) forms the entire state space
const states_2d = [(x, v) for x in x1 for v in x2]
const nstates = length(states_2d)

# Action discretization
const actions = collect(LinRange(u_min, u_max, num_points_action))
const nactions = length(actions)

println("Number of states  = $nstates")
println("Number of actions = $nactions")

##########################################################
# 2) Build the transition probabilities T[s, a, s_next]
##########################################################
# T is a 3D array of size (nstates × nactions × nstates),
# T[s, a, s_next] = Probability of going to s_next from s under action a.
# The task is to keep the state trajectory inside the box K = [0, 4] × [−3, 3],
# thus the target is taken to be its complement T = KC.

function is_safe(x::Float64, v::Float64)
    # "Safe" region: x in [0,4], v in [-3,3]
    return (0.0 <= x <= 4.0) && (-3.0 <= v <= 3.0)
end

# Continuous (noisy) dynamics for double-integrator
function dynamics_rand(x::Float64, v::Float64, u::Float64)    
    if !is_safe(x, v)
        return (x, v)   # no movement if outside the safe region
    else
        # make sure d is between d min and d max 
        d = rand(Normal(0, sigma))
    
        # Ensure d is between d_min and d_max
        while d < d_min || d > d_max
            d = rand(Normal(0, sigma))
        end
        dt=0.1
        x1_next = x + v*dt
        x2_next = v + (u+d)*dt
        return (x1_next, x2_next)
    end
end

# Build a KD-tree for snapping continuous next-states to the nearest discrete state
states_matrix = hcat([collect(s) for s in states_2d]...)  # shape: 2 x nstates
tree = KDTree(states_matrix)

println("\nBuilding transition probabilities T[s, a, s_next] ...")

# Initialize T array: T[s, a, s_next]
T = zeros(Float64, nstates, nactions, nstates)

# (Compute Transition Matrix)
for is in 1:nstates
    s = states_2d[is]   # e.g., s = (x, v)
    for a in 1:nactions
        d = rand(Normal(0, sigma))
        for i in 1:nsamples
            xn = s[1]
            vn = s[2]
            j = 1
            while true	
                (xn, vn) = dynamics_rand(s[1], s[2], actions[a])
                # For knn, pass a 2-element Vector, not a Tuple
                idxs, dists = knn(tree, [xn, vn], 1)
                if first(idxs) != is || j > 1
                    T[is, a, first(idxs)] += 1.0 / nsamples
                    break
                end
                j += 1
            end
        end
    end
end

@assert sum(T[1,1,:]) ≈ 1.0  # Checks row sums to ~1
@assert all(T .>= 0.0) "Matrix T contains negative values!"

for action in 1:num_points_action
    @assert all(i -> sum(T[i, action, :]) ≈ 1.0, axes(T, 1)) "Not all row sums are approximately 1!"
end

# Remove very small transition probabilities
T .= ifelse.(abs.(T) .< threshold_for_transit, 0.0, T)

println("Done building T.\n")

##########################################################
# 3) Reward function
##########################################################
# +1 if in safe region, else 0
function reward(s::Tuple{Float64, Float64})
    return is_safe(s[1], s[2]) ? 1.0 : 0.0
end

# Initialize reward vector for all states
r = zeros(nstates)
for s in 1:nstates
    r[s] = reward(states_2d[s])
end

##########################################################
# 4) Solve the linear program (primal formulation)
##########################################################
function solve_primal_case()   
    # Setup model
    model = Model(Gurobi.Optimizer)
    
    # Define variables: g[s] and h[s] for each state
    @variable(model, g[1:nstates] >= 0)
    @variable(model, h[1:nstates] >= 0)
    
    # Define alpha distribution (uniform for this example)
    alpha_dist = ones(nstates) ./ nstates
    
    # Objective: minimize sum_{s} alpha_s * g_s
    @objective(model, Min, sum(alpha_dist[s] * g[s] for s in 1:nstates))
    
    # Constraint 1: g_s ≥ sum_{s'} g_s' * P(s'|s,a) for all s in S, a in A
    for s in 1:nstates
        for a in 1:nactions
            @constraint(model, g[s] >= sum(g[s_prime] * T[s, a, s_prime] for s_prime in 1:nstates))
        end
    end
    
    # Constraint 2: g_s + h_s ≥ r(s) + sum_{s'} h_s' * P(s'|s,a) for all s in S, a in A
    for s in 1:nstates
        for a in 1:nactions
            @constraint(model, g[s] + h[s] >= r[s] + sum(h[s_prime] * T[s, a, s_prime] for s_prime in 1:nstates))
        end
    end

    # Solve the model
    optimize!(model)
    
    # Check solution status
    stat = termination_status(model)
    println("Solver status: $stat")
    
    if stat == MOI.OPTIMAL
        println("Optimal objective value = ", objective_value(model))
        
        # Retrieve the optimal g and h values
        g_opt = value.(g)
        h_opt = value.(h)
        
        # Reshape for visualization
        g_map = reshape(g_opt, num_points_state, num_points_state)
        h_map = reshape(h_opt, num_points_state, num_points_state)
        
        # Create grid for visualization
        X = [x1[i] for i in 1:num_points_state, j in 1:num_points_state]
        Y = [x2[j] for i in 1:num_points_state, j in 1:num_points_state]
        
        return objective_value(model), g_opt, h_opt, g_map, h_map, X, Y
    else
        println("No optimal solution found. Status = ", stat)
        return nothing, nothing, nothing, nothing, nothing, nothing, nothing
    end
end

##########################################################
# 5) Save results
##########################################################
function main_primal()
    # Solve the optimization problem
    println("*************************************************")
    println("Solving primal case with sigma = $sigma")
    println("*************************************************")
    
    # Call solve_primal_case and capture the returned values
    objective, g_opt, h_opt, g_map, h_map, X, Y = solve_primal_case()
    
    if isnothing(objective)
        println("Failed to solve the optimization problem.")
        return
    end
    
    # Create a timestamped folder for results
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    results_dir = joinpath(@__DIR__, "results_primal", string(sigma))
    mkpath(results_dir)
    
    println("Saving results to: $results_dir")
    
    # Save optimization results
    open(joinpath(results_dir, "optimization_results.txt"), "w") do f
        println(f, "Optimization run on: $timestamp")
        println(f, "Sigma: $sigma")
        println(f, "Objective value (min average reward): $objective")
    end
    
    # Save g and h vectors
    writedlm(joinpath(results_dir, "g_vector.csv"), g_opt, ',')
    writedlm(joinpath(results_dir, "h_vector.csv"), h_opt, ',')
    
    # Save g and h maps
    writedlm(joinpath(results_dir, "G.csv"), g_map, ',')
    writedlm(joinpath(results_dir, "H.csv"), h_map, ',')
    
    # Save grid data
    writedlm(joinpath(results_dir, "X_grid.csv"), X, ',')
    writedlm(joinpath(results_dir, "Y_grid.csv"), Y, ',')
    
    # Extract optimal policy (for each state, find action that maximizes the value)
    optimal_policy = zeros(Int, nstates)
    for s in 1:nstates
        max_val = -Inf
        best_a = 1
        
        for a in 1:nactions
            val = r[s] + sum(h_opt[s_prime] * T[s, a, s_prime] for s_prime in 1:nstates)
            if val > max_val
                max_val = val
                best_a = a
            end
        end
        
        optimal_policy[s] = best_a
    end
    
    # Save optimal policy
    policy_map = reshape(optimal_policy, num_points_state, num_points_state)
    writedlm(joinpath(results_dir, "optimal_policy.csv"), policy_map, ',')
    
    println("Results saved successfully to: $results_dir")
end

##########################################################
# 6) Run
##########################################################
main_primal()
