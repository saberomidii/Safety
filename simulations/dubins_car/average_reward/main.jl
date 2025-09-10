##########################################################
# Average Reward and Linear Programming for Inverted Pendulum Systems 
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

using MosekTools       
import MathOptInterface as MOI
using Mosek

using SparseArrays
using LinearAlgebra


##########################################################
# 1) Setup: parameters, discretization, etc.
##########################################################
# Optional: set a random seed if desired
Random.seed!(10)

# Single noise standard deviation
const sigma = 1.0
const threshold_for_transit = 0.001

const x_1_min = -2.0
const x_1_max = 2.0

const x_2_min = -2.0
const x_2_max =  2.0

const x_3_min = 0
const x_3_max = 2π

const u_min = -1.0
const u_max =  1.0

const d_min=-1
const d_max= 1

const num_points_action =21

const num_points_state_1=101
const num_points_state_2=101
const num_points_state_3=18



# Number of random samples used when constructing transitions
const nsamples = 100

# Number of states in each dimension
x1 = collect(LinRange(x_1_min, x_1_max, num_points_state_1))  # Possible x values
x2 = collect(LinRange(x_2_min, x_2_max, num_points_state_2))  # Possible v values
x3 = collect(LinRange(x_3_min, x_3_max, num_points_state_3))  # Possible v values

# Cartesian product of (x, v) forms the entire state space
const states_3d = [(x, v, theta) for x in x1 for v in x2 for theta in x3]
const nstates = length(states_3d)

# Action discretization
const actions = collect(LinRange(u_min, u_max, num_points_action))
const nactions = length(actions)

println("Number of states  = $nstates")
println("Number of actions = $nactions")

##########################################################
# 2) Build the transition probabilities T[s, a, s_next]
# ##########################################################

function is_safe(x::Float64, v::Float64)
    # Center of both ellipses
    x_c, v_c = 0.0, 0.0

    # Larger ellipse (outer boundary)
    a_outer, b_outer = 0.4, 0.9

    # Smaller ellipse (inner hole)
    a_inner, b_inner = 0.2, 0.70

    # Check if point is inside the larger ellipse
    inside_outer = ((x - x_c)^2) / a_outer^2 + ((v - v_c)^2) / b_outer^2 <= 1

    # Check if point is inside the smaller ellipse
    inside_inner = ((x - x_c)^2) / a_inner^2 + ((v - v_c)^2) / b_inner^2 <= 1

    # Safe set is the ring: inside the outer, but outside the inner
    return inside_outer && !inside_inner
end

# Continuous (noisy) dynamics for Inverted Pendulum
function dynamics_rand(x::Float64, y::Float64, theta::Float64, u::Float64)    
    if !is_safe(x, y)
        return (x, y, theta)   # no movement if outside the safe region
    else
        # make sure d is between d min and d max 
        d = rand(Normal(0, sigma))
    
        # Ensure d is between d_min and d_max
        while d < d_min || d > d_max
            d = rand(Normal(0, sigma))
        end
        dt=0.1
        V=0.1       #Constant Speed 
        x1_next = x + V*cos(theta)*dt
        x2_next = y + V*cos(theta)*dt
        x3_next = theta + (u+d)*dt
        return (x1_next, x2_next, x3_next)
    end
end

# Build a KD-tree for snapping continuous next-states to the nearest discrete state
states_matrix = hcat([collect(s) for s in states_3d]...)  # shape: 3 x nstates
tree = KDTree(states_matrix)

println("\nBuilding transition probabilities T[s, a, s_next] ...")

# Initialize T array: T[s, a, s_next]
T = Vector{SparseMatrixCSC{Float64,Int64}}(undef,nstates)

for index_state in 1:nstates
	T[index_state] = spzeros(nactions,nstates)
end

# (Compute Transition Matrix)
@time begin
 for is in 1:nstates
    s = states_3d[is]   # e.g., s = (x, y, theta)
    for a in 1:nactions
        d = rand(Normal(0, sigma))
        for i in 1:nsamples
            xn = s[1]
            vn = s[2]
            thetan=s[3]
            j = 1
            while true	
                (xn, vn, thetan) = dynamics_rand(s[1], s[2],s[3], actions[a])
                # For knn, pass a 2-element Vector, not a Tuple
                idxs, dists = knn(tree, [xn, vn, thetan], 1)
                if first(idxs) != is || j > 1
		   T[is][ a, first(idxs)] += 1.0 / nsamples
                   break
                end
                j += 1
            end
        end
    end
 end
end

num_zeros = count(x -> x==0,T)

println("Negative Values checking") 
@assert all(t->all(t.>=0.0),T) "Negative Value!"

println("Rows sums are approximately 1!")

for action in 1:num_points_action
	@assert all(i -> sum(T[i][action,:])≈1.0, axes(T,1)) "Not All row sums are approximately 1!"
end

println("Removing small values")
for s in 1:nstates
	T[s]= dropzeros!(map(x->abs(x)<threshold_for_transit ? 0.0 : x, T[s]))
end

println("Done Building T.\n")

##########################################################
# 3) Reward function
##########################################################
# +1 if in safe region, else 0
function reward(s::Tuple{Float64, Float64, Float64})
    return is_safe(s[1], s[2]) ? 1.0 : 0.0
end

# Initialize reward vector for all states
r = zeros(nstates)
for s in 1:nstates
    r[s] = reward(states_3d[s])
end

##########################################################
# 4) Solve the linear program (primal formulation)
##########################################################
function solve_primal_case()   
    # Setup model
    model = Model(Mosek.Optimizer)

    # --- CORRECTED TOLERANCE SETTINGS FOR A LINEAR PROGRAM ---
    # Set the tolerance for the interior-point optimizer's relative gap.
    set_optimizer_attribute(model, "MSK_DPAR_INTPNT_CO_TOL_REL_GAP", 1e-4)

    # You can also set primal and dual feasibility tolerances if needed.
    set_optimizer_attribute(model, "MSK_DPAR_INTPNT_CO_TOL_PFEAS", 1e-4)
    set_optimizer_attribute(model, "MSK_DPAR_INTPNT_CO_TOL_DFEAS", 1e-4)
    
    # Setting this to 0 (MSK_OFF) disables the basis identification procedure.
    set_optimizer_attribute(model, "MSK_IPAR_INTPNT_BASIS", 0)
    
    set_time_limit_sec(model, 1800.0)

    # Define variables: g[s] and h[s] for each state
    @time begin
    @variable(model, g[1:nstates] >= 0)
    @variable(model, h[1:nstates] >= 0)
    end

    # Define alpha distribution (uniform for this example)
    alpha_dist = ones(nstates) ./ nstates
    
    # Objective: minimize sum_{s} alpha_s * g_s
    @time begin
    @objective(model, Min, sum(alpha_dist[s] * g[s] for s in 1:nstates))
    end


    # Precompute the indices of nonzero transitions for each (s,a)
    nonzero_indices = Vector{Vector{Int}}(undef, nstates * nactions)

    for s in 1:nstates
        for a in 1:nactions
            idx = (s - 1) * nactions + a
            # Find all s_prime where T[s,a,s_prime] > threshold
	    nz_sprimes = findall(p -> p > threshold_for_transit, T[s][a,:])
            nonzero_indices[idx] = nz_sprimes
        end
    end

    @time begin
    for s in 1:nstates
        for a in 1:nactions
            idx = (s - 1) * nactions + a
            s_primes = nonzero_indices[idx]
            if !isempty(s_primes)
		    @constraint(model, g[s] >= sum(g[s_prime] * T[s][ a, s_prime] for s_prime in s_primes))
            end
        end
    end
    end

    @time begin
    for s in 1:nstates
        for a in 1:nactions
            idx = (s - 1) * nactions + a
            s_primes = nonzero_indices[idx]
            if !isempty(s_primes)
		    @constraint(model, g[s] + h[s] >= r[s] + sum(h[s_prime] * T[s][a, s_prime] for s_prime in s_primes))
            end
        end
    end
    end

    # Solve the model with error handling
    try
        # print_model_summary(model)
        optimize!(model)
        
        # Check solution status
        stat = termination_status(model)
        println("Solver status: $stat")
        
        if stat == MOI.OPTIMAL || stat == MOI.ALMOST_OPTIMAL || stat == MOI.ITERATION_LIMIT || stat == MOI.TIME_LIMIT
            println("Solution found with status: $stat")
            println("Objective value = ", objective_value(model))
            
            # Retrieve the g and h values (optimal or best found)
            g_opt = value.(g)
            h_opt = value.(h)
            
            # Check for NaN or Inf values and replace them
            g_opt = replace(g_opt, NaN => 0.0, Inf => 1.0, -Inf => 0.0)
            h_opt = replace(h_opt, NaN => 0.0, Inf => 1.0, -Inf => 0.0)
            
            # Round to two decimal places
            g_opt = round.(g_opt, digits=2)
            h_opt = round.(h_opt, digits=2)
            
            # Reshape for visualization
            g_map = reshape(g_opt, num_points_state_1, num_points_state_2, num_points_state_3)
            h_map = reshape(h_opt, num_points_state_1, num_points_state_2, num_points_state_3)
            
            # Create grid for visualization
            X = [x1[i] for i in 1:num_points_state_1, j in 1:num_points_state_2]
            Y = [x2[j] for i in 1:num_points_state_1, j in 1:num_points_state_2]
            
            return objective_value(model), g_opt, h_opt, g_map, h_map, X, Y
        else
            println("No optimal solution found. Status = ", stat)
            return nothing, nothing, nothing, nothing, nothing, nothing, nothing
        end
    catch e
        println("Error during optimization: ", e)
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
    # results_dir = joinpath(@__DIR__, "results_primal", string(sigma))
    results_dir = joinpath(@__DIR__, "results_primal_$(sigma)")

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
    @time begin
    for s in 1:nstates
        max_val = -Inf
        best_a = 1
        
        for a in 1:nactions
		val = r[s] + sum(h_opt[s_prime] * T[s][a, s_prime] for s_prime in 1:nstates)
            if val > max_val
                max_val = val
                best_a = a
            end
        end
        
        optimal_policy[s] = best_a
    end
    end 

    # Save optimal policy
    policy_map = reshape(optimal_policy, num_points_state_1, num_points_state_2, num_points_state_3)
    writedlm(joinpath(results_dir, "optimal_policy.csv"), policy_map, ',')
    
    println("Results saved successfully to: $results_dir")
end

##########################################################
# 6) Run
##########################################################
main_primal()

 


