##########################################################
# Average Reward and Linear Programming for Double-Integrator Systems 
# A Minimum Discounted Reward Hamilton–Jacobi Formulation for Computing Reachable Sets
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
# Set Gurobi WLS license credentials - add these at the to
# # At the very beginning of your file, before any other code executes:
#println("Setting Gurobi license credentials...")
ENV["GRB_WLSACCESSID"] = "52eb20bf-115c-42d3-931f-47561460611c"
ENV["GRB_WLSSECRET"] = "587b8f93-6d53-43c9-af49-6b96ac589004"
ENV["GRB_LICENSEID"] = "2611020"
println("License credentials set")
# Optional: set a random seed if desired
Random.seed!(50)

# Single noise standard deviation
<<<<<<< HEAD
const sigma = 2.0
=======
const sigma = 0.35
>>>>>>> 223909fd30b74269165535b434ed30d5f0e43367

const treshold_for_transit=0.001
 
const x_1_min=- 1.0
const x_1_max=  5.0
const x_2_min= -5.0
const x_2_max=  5.0

const u_min= -2.0
const u_max=  2.0


<<<<<<< HEAD
const num_points_action = 5
=======
const num_points_action = 11
>>>>>>> 223909fd30b74269165535b434ed30d5f0e43367

const num_points_state = 55

const num_points_state = 55

# Number of random samples used when constructing transitions
const nsamples = 100

# Number of states in each dimension
x1 = collect(LinRange(x_1_min, x_2_max, num_points_state))  # Possible theta values
x2 = collect(LinRange(x_2_min, x_2_max, num_points_state))  # Possible theta_dot values

# Cartesian product of (x, v) forms the entire state space
const states_2d  = [(x, v) for x in x1 for v in x2]
const nstates    = length(states_2d)


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
# thus the target is taken to be its complement T = KC . For ease of 
# exposition we define the safe set Ω(T ) := R(T )C .

function is_safe(x::Float64, v::Float64)
    # Example "safe" region: x in [-1,1], v in [-1,1]
    return (0.0<= x <= 4.0) && (-3.0<= v <= 3.0)
end

# Continuous (noisy) dynamics:
## inverted pendulumn 
function dynamics_rand(x::Float64, v::Float64, u::Float64,d)    
    if !is_safe(x, v)
        return (x, v)   # no movement if outside the safe region
    else
        # d=rand(Normal(0,sigma))
        dt=0.1
	d=rand(Normal(0,sigma))
        x1_next = x + v*dt
        x2_next = v + (u+d)*dt
        return (x1_next, x2_next)
    end
end

# Build a KD-tree for snapping continuous next-states to the nearest discrete state.
states_matrix = hcat([collect(s) for s in states_2d]...)  # shape: 2 x nstates
tree = KDTree(states_matrix)

println("\nBuilding transition probabilities T[s, a, s_next] ...")

# Initialize T array: T[s, a, s_next]
T = zeros(Float64, nstates, nactions, nstates)
disturbance_list=zeros(nstates, nactions, nsamples)


d = [rand(Normal(0, sigma)) for _ in 1:nsamples]

# (Compute Transition Matrix)
for is in 1:nstates
    s = states_2d[is]   # e.g., s = (x, v)
    for a in 1:nactions
        for i in 1:nsamples
	    xn=s[1]
	    vn=s[2]
	    j=1
	    while true	
            (xn, vn) = dynamics_rand(s[1], s[2], actions[a],d[i])
            # (xn, vn) = dynamics_rand(s[1], s[2], actions[a])
            # For knn, pass a 2-element Vector, not a Tuple
	    idxs, dists = knn(tree, [xn, vn], 1)
	             if first(idxs) != is || j > 1 
		                  T[is, a, first(idxs)] += 1.0 / nsamples
		                  break
	            end   # Create meshgrid-like structure end
	    j+=1
            end
        end
    end
end

@assert sum(T[1,1,:]) ≈ 1.0  # Checks row sums to ~1
@assert all(T .>= 0.0) "Matrix T contains negative values!"

for action in 1:num_points_action
@assert all(i -> sum(T[i, action, :]) ≈ 1.0, axes(T, 1)) "Not all row sums are approximately 1!"
end

T .= ifelse.(abs.(T) .< treshold_for_transit, 0.0, T)

println("The maximum disturbance value is ",maximum(disturbance_list))

println("The minumum disturbance value is ",minimum(disturbance_list))




println("Done building T.\n")


##########################################################
# 3) Reward function
##########################################################
# +1 if in safe region, else 0
function reward(s::Tuple{Float64, Float64}, a::Float64)
    return is_safe(s[1], s[2]) ? 1.0 : 0.0
end

##########################################################
# 4) Solve the linear program (flow constraints)
##########################################################
function solve_case()   
    # Set Gurobi license credentials right before creating the model
   # ENV["GRB_WLSACCESSID"] = "52eb20bf-115c-42d3-931f-47561460611c"
   # ENV["GRB_WLSSECRET"] = "587b8f93-6d53-43c9-af49-6b96ac589004"
   # ENV["GRB_LICENSEID"] = "2611020"

    # Setup model
    model = Model(Gurobi.Optimizer)

    # Define z[s, a] and y[s, a]
    @variable(model, z[1:nstates, 1:nactions] >= 0)
    @variable(model, y[1:nstates, 1:nactions] >= 0)

    # c1: total measure of z is 1
    @constraint(model, c1, sum(z[s,a] for s in 1:nstates, a in 1:nactions) == 1)

    # c2: flow constraints
    @constraint(model, c2[j in 1:nstates],
        sum(z[s,a] * T[s,a,j] for s in 1:nstates, a in 1:nactions)
        == sum(z[j,a2]       for a2 in 1:nactions)
    )

    # c3: alpha distribution
    # For demonstration, alpha = 1 / nstates (uniform)
    alpha_dist = 1.0 / nstates
    @constraint(model, c3[j in 1:nstates],
        sum(z[j,a] for a in 1:nactions)
      + sum(y[j,a] for a in 1:nactions)
      - sum(y[s,a] * T[s,a,j] for s in 1:nstates, a in 1:nactions)
        == alpha_dist
    )

    # Objective: maximize sum_{s,a} z[s,a] * reward(s,a)
    @objective(model, Max, sum(z[s,a]*reward(states_2d[s], actions[a])
                               for s in 1:nstates, a in 1:nactions))

    optimize!(model)
    stat = termination_status(model)
    println("Solver status: $stat")
    if stat == MOI.OPTIMAL
        println("Optimal objective value = ", objective_value(model))
    else
        println("No optimal solution found. status = ", stat)
    end

    # Retrieve primal solution: z_sol[s] = sum_{a} z[s,a]
    z_sol = zeros(nstates)
    for s in 1:nstates
        z_sol[s] = sum(value(z[s,a]) for a in 1:nactions)
    end

    # Build Nx×Ny matrix for occupation measure z(s).
    Nx = length(x1)
    Ny = length(x2)
    Z_occup = zeros(Nx, Ny)

    # The indexing in states_2d is row-major: 
    # states_2d[(i-1)*Ny + j] = (positions[i], velocities[j])
    for i in 1:Nx
        for j in 1:Ny
            sidx = (i - 1)*Ny + j
            Z_occup[i,j] = z_sol[sidx]
        end
    end

    # Prepare meshgrid arrays (X,Y) for plotting
    X = [x1[i]  for j in 1:Ny, i in 1:Nx]  # Nx×Ny
    Y = [x2[j] for j in 1:Ny, i in 1:Nx]

    # Retrieve duals for c2, c3
    dual_c2_vals = [dual(c2[j]) for j in 1:nstates]
    dual_c3_vals = [dual(c3[j]) for j in 1:nstates]

    # Reshape them into Nx×Ny
    G_map = zeros(Nx, Ny)
    H_map = zeros(Nx, Ny)
    for i in 1:Nx
        for j in 1:Ny
            sidx = (i - 1)*Ny + j
            G_map[i,j] = -dual_c2_vals[sidx]
            H_map[i,j] = -dual_c3_vals[sidx]
        end
    end

    return objective_value(model), X, Y, Z_occup, G_map, H_map
end

##########################################################
# 5) Plot in 3D (z, g, h) side-by-side (No Colorbars)
#    Then save results to a new directory "results/"
##########################################################
function main_3D()
    # Solve the optimization problem
    println("*************************************************")
    println("Solving case with sigma = $sigma")
    println("*************************************************")
    
    # Call solve_case() and capture all the returned values
   
    objective, X, Y, Z_occup, G_map, H_map = solve_case()
    
    # Extract any other values you need for results
    
    # Create a timestamped folder for results
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    results_dir = joinpath(@__DIR__, "results_dual", string(sigma))
    mkpath(results_dir)
    
    println("Saving results to: $results_dir")
    
    # Save optimization results
    open(joinpath(results_dir, "optimization_results.txt"), "w") do f
        println(f, "Optimization run on: $timestamp")
        println(f, "Sigma: $sigma")
        println(f, "Average reward: $objective")
        # Add other parameters you want to save
    end
    
    # Save occupancy measures data
    if !isnothing(Z_occup)
        writedlm(joinpath(results_dir, "occup_measures_Z_occup.csv"), Z_occup, ',')
    end
    
    # Save value function data if available
    if !isnothing(G_map)
        writedlm(joinpath(results_dir, "G_map.csv"), G_map, ',')
    end
    
    if !isnothing(H_map)
	    writedlm(joinpath(results_dir, "H_map.csv"), H_map, ',')
    end

    # Construct grid data for plotting elsewhere if needed
    #X, Y = construct_grid()  # You'll need to implement this function
    writedlm(joinpath(results_dir, "X_grid.csv"), X, ',')
    writedlm(joinpath(results_dir, "Y_grid.csv"), Y, ',')
    
    println("Results saved successfully to: $results_dir")
end
# Add this function to your code to fix the error:

##########################################################
# Helper function to construct the grid for saving results
##########################################################
function construct_grid()
    # Use the same grid parameters as in your original problem
    nx = num_points_state  # Use the same values as in your state space discretization
    ny = num_points_state
    
    # These ranges should match your state space

    xmin, xmax = x_1_min, x_1_max
    ymin, ymax = -x_2_min, x_2_max
    
    # Create a grid
    x = range(xmin, xmax, length=nx)
    y = range(ymin, ymax, length=ny)
    

    X = [i for i in x, j in 1:ny]
    Y = [j for i in 1:nx, j in y]
    
    return X, Y
end
##########################################################
# 6) Run
##########################################################
main_3D()
