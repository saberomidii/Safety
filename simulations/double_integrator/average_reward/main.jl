##########################################################
# Average Reward and Linear Programming for Double-Integrator Systems 
# Primal Formulation for Computing Safe Sets
##########################################################
using Dates
using DelimitedFiles  # For saving data to CSV files
using Distributions
using StatsFuns
using Statistics
using Random
using NearestNeighbors
using JuMP
using FileIO   # for saving results if needed
using Base: mkpath


using MosekTools
using Mosek
import MathOptInterface as MIO


using SparseArrays
using LinearAlgebra 

##########################################################
# 1) Setup: parameters, discretization, etc.
##########################################################

Random.seed!(10) # to reproduce the results. 

# Single noise standard deviation
const sigma = 1.0

const threshold_for_transit = 0.001

const x_1_min = -1.0
const x_1_max = 5.0
const x_2_min = -5.0
const x_2_max = 5.0

const u_min = -2.0
const u_max = 2.0

const d_min= -1.0
const d_max=  1.0 

const num_points_action = 41
  
const num_points_state_1 = 161
const num_points_state_2 = 161

# Number of random samples used when constructing transitions
const nsamples = 100

# Number of states in each dimension
x1 = collect(LinRange(x_1_min, x_1_max, num_points_state_1))  # Possible x values
x2 = collect(LinRange(x_2_min, x_2_max, num_points_state_2))  # Possible v values

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

println("\nBuilding transition probabilities T[s][a, s_next] ...")

# Initialize T array: T[s][a, s_next]
T= Vector{SparseMatrixCSC{Float64,Int64}}(undef,nstates)

for index_state in 1:nstates
	T[index_state] = spzeros(nactions,nstates)
end

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
			T[is][a, first(idxs)] += 1.0 / nsamples
                    break
                end
                j += 1
            end
        end
    end
end

println("Negative values check")
@assert all(t->all(t.>=0.0),T) "Negative Value!"

println("Rows sums are approximately 1!")
for action in 1:num_points_action
	@assert all(i -> sum(T[i][action, :]) ≈ 1.0, axes(T, 1)) "Not All row sums are approximately1!"
end

println("Removing small values")
for s in 1:nstates
	T[s]= dropzeros!(map(x->abs(x)<threshold_for_transit ? 0.0 : x, T[s]))
end


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
    model = Model(Mosek.Optimizer)
   
    # --- CORRECTED TOLERANCE SETTINGS FOR A LINEAR PROGRAM ---
    # Set the tolerance for the interior-point optimizer's relative gap.
    set_optimizer_attribute(model, "MSK_DPAR_INTPNT_CO_TOL_REL_GAP", 1e-5)

    # You can also set primal and dual feasibility tolerances if needed.
    set_optimizer_attribute(model, "MSK_DPAR_INTPNT_CO_TOL_PFEAS", 1e-5)
    set_optimizer_attribute(model, "MSK_DPAR_INTPNT_CO_TOL_DFEAS", 1e-5)
    
    # Setting this to 0 (MSK_OFF) disables the basis identification procedure.
    set_optimizer_attribute(model, "MSK_IPAR_INTPNT_BASIS", 0)
    
    set_time_limit_sec(model, 15*60)

    
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

    # Precompute nonzero for s and a 
    nonzero_indices = Vector{Vector{Int64}}(undef,nstates*nactions)

    for s in 1:nstates
	    for a in 1:nactions
		idx = (s-1)*nactions + a

		nz_sprimes =findall(p -> p> threshold_for_transit,T[s][a,:])
		nonzero_indices[idx] = nz_sprimes
	    end 
    end

    # Constraint 1: g_s ≥ sum_{s'} g_s' * P(s'|s,a) for all s in S, a in A
    @time begin 
    		for s in 1:nstates
        		for a in 1:nactions
			    idx = (s-1)*nactions + a
			    s_primes= nonzero_indices[idx]

            			if !isempty(s_primes)
					@constraint(model, g[s] >= sum(g[s_prime]*T[s][a, s_prime] for s_prime in s_primes))
            			end
        		end
   		 end
	 end
    # Constraint 2: g_s + h_s ≥ r(s) + sum_{s'} h_s' * P(s'|s,a) for all s in S, a in A
    @time begin 
	       for s in 1:nstates
        		for a in 1:nactions

			       idx = (s-1)*nactions + a 
			       s_primes = nonzero_indices[idx]
			       if !isempty(s_primes)
				       @constraint(model, g[s]+h[s] >= r[s] + sum(h[s_prime]*T[s][a,s_prime] for s_prime in s_primes))
			       end
			end
	   	end 
	   end

    # Solve the model with error handling
    try
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
            g_map = reshape(g_opt, num_points_state_1, num_points_state_2)
            h_map = reshape(h_opt, num_points_state_1, num_points_state_2)
            
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
		val = r[s] + sum(h_opt[s_prime] * T[s][a, s_prime] for s_prime in 1:nstates)
            if val > max_val
                max_val = val
                best_a = a
            end
        end
        
        optimal_policy[s] = best_a
    end
    
    # Save optimal policy
    policy_map = reshape(optimal_policy, num_points_state_1, num_points_state_2)
    writedlm(joinpath(results_dir, "optimal_policy.csv"), policy_map, ',')
    
    println("Results saved successfully to: $results_dir")
end

##########################################################
# 6) Run
##########################################################
main_primal()
