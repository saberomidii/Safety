##########################################################
# Average Reward and Linear Programming for Double-Integrator Systems
# GPU-Optimized Version with Visualization
##########################################################
using Distributions
using StatsFuns
using Statistics
using Random
using NearestNeighbors
using JuMP
import Gurobi
using GLMakie
using FileIO   # for GLMakie.save if needed
using Base: mkpath
using CUDA

# Check if CUDA is available
const has_cuda = CUDA.functional()
if has_cuda
    println("CUDA is available! GPU acceleration enabled.")
else
    println("CUDA is not available. Falling back to CPU implementation.")
end

##########################################################
# 1) Setup: parameters, discretization, etc.
##########################################################
# At the beginning of your file, set the Gurobi license credentials
println("Setting Gurobi license credentials...")
ENV["GRB_WLSACCESSID"] = "52eb20bf-115c-42d3-931f-47561460611c"
ENV["GRB_WLSSECRET"] = "587b8f93-6d53-43c9-af49-6b96ac589004"
ENV["GRB_LICENSEID"] = "2611020"
println("License credentials set")

# Optional: set a random seed for reproducibility
# Random.seed!(42)

# Parameters
const sigma = 0.2
const trashhold_for_transit = 0.001

const theta_min = -0.35
const theta_max = 0.35
const theta_dot_min = -0.65
const theta_dot_max = 0.65

const u_min = -3.0
const u_max = 3.0

const num_points_action = 10
const num_points_state = 51

# Number of random samples used when constructing transitions
const nsamples = 1000

# Number of states in each dimension
theta = collect(LinRange(theta_min, theta_max, num_points_state))
theta_dot = collect(LinRange(theta_dot_min, theta_dot_max, num_points_state))

# Cartesian product of (theta, theta_dot) forms the entire state space
const states_2d = [(x, v) for x in theta for v in theta_dot]
const nstates = length(states_2d)

# Action discretization
const actions = collect(LinRange(u_min, u_max, num_points_action))
const nactions = length(actions)

println("Number of states  = $nstates")
println("Number of actions = $nactions")

##########################################################
# 2) Build the transition probabilities T[s, a, s_next]
##########################################################
# Safe region check - used in dynamics
function is_safe(x::Float64, v::Float64)
    return (-0.30 <= x <= 0.30) && (-0.60 <= v <= 0.60)
end

# Continuous (noisy) dynamics: inverted pendulum
function dynamics_rand(theta::Float64, theta_dot::Float64, u::Float64)
    d = rand(Normal(0, sigma))
    # d = sin(randn())
    if !is_safe(theta, theta_dot)
        return (theta, theta_dot, d)   # no movement if outside the safe region
    else
    # Physical parameters of the inverted pendulum
    g = 9.81   # Gravity (m/s²)
    m = 2.0    # Mass (kg)
    l = 1.0    # Length (m)
    b = 0.1    # Damping coefficient (kg m²/s) (viscous friction)
    I = m * l^2  # Moment of inertia (kg·m²)
    dt = 0.1  # Time step (s)
    
        theta_next = theta + theta_dot*dt
        theta_dot_next = theta_dot + ((g/l*sin(theta) + (1/m*l^2)*u) + d)*dt
        # theta_dot_next = theta_dot + ( (g / l) * sin(theta) - (b / I) * theta_dot + (1 / I) * u + d ) * dt

        return (theta_next, theta_dot_next, d)
    end
end

# function dynamics_rand(theta::Float64, theta_dot::Float64, u::Float64)
#     d = rand(Normal(0, sigma))  # Stochastic disturbance (Gaussian noise)

#     if !is_safe(theta, theta_dot)
#         return (theta, theta_dot, d)  # No movement if outside the safe region
#     else
#         # Equations of motion for inverted pendulum
#         theta_next = theta + theta_dot * dt
#         theta_dot_next = theta_dot + ( (g / l) * sin(theta) - (b / I) * theta_dot + (1 / I) * u + d ) * dt
        
#         return (theta_next, theta_dot_next, d)
#     end
# end



# GPU-optimized function to compute transitions in batches if CUDA is available
function build_transition_matrix()
    println("\nBuilding transition probabilities T[s, a, s_next] ...")
    
    # Initialize T array: T[s, a, s_next]
    T = zeros(Float64, nstates, nactions, nstates)
    disturbance_list = zeros(nstates, nactions, nsamples)
    
    # Build a KD-tree for snapping continuous next-states to the nearest discrete state
    states_matrix = hcat([collect(s) for s in states_2d]...)  # shape: 2 x nstates
    tree = KDTree(states_matrix)
    
    # Use parallel processing if available
    # if has_cuda && nstates * nactions * nsamples > 10^6
        # This is a placeholder for GPU implementation
        # In a real implementation, you would use CUDA kernels to compute transitions in parallel
        # For now we'll use CPU implementation with threading
        Threads.@threads for is in 1:nstates
            s = states_2d[is]
            for a in 1:nactions
                for i in 1:nsamples
                    xn=s[1]
                    vn=s[2]
                    j=1
                    while true 
                        (xn, vn, disturbance) = dynamics_rand(xn, vn, actions[a])
                        # disturbance_list[is, a, i] = disturbance
                        # Find nearest state
                        idxs, dists = knn(tree, [xn, vn], 1)
                        if first(idxs) != is || j >10
                            T[is, a, first(idxs)] += 1.0 / nsamples
                            break
                        end
                        j+=1
                    end
                    
                end
            end
        end
    # else
        # # CPU implementation with threading
        # Threads.@threads for is in 1:nstates
        #     s = states_2d[is]
        #     for a in 1:nactions
        #         for i in 1:nsamples
        #             (xn, vn, disturbance) = dynamics_rand(s[1], s[2], actions[a])
        #             disturbance_list[is, a, i] = disturbance
        #             # Find nearest state
        #             idxs, dists = knn(tree, [xn, vn], 1)
        #             T[is, a, first(idxs)] += 1.0 / nsamples
                # end
            # end
        # end
    # end
    
    # Clean up small values for memory efficiency
    T .= ifelse.(abs.(T) .< trashhold_for_transit, 0.0, T)
    
    println("The maximum disturbance value is ", maximum(disturbance_list))
    println("The minimum disturbance value is ", minimum(disturbance_list))
    println("Done building T.\n")
    
    # Verify correctness
    @assert sum(T[1, 1, :]) ≈ 1.0 "First row doesn't sum to 1"
    @assert all(T .>= 0.0) "Matrix T contains negative values!"
    
    for action in 1:num_points_action
        @assert all(i -> sum(T[i, action, :]) ≈ 1.0, axes(T, 1)) "Not all row sums are approximately 1!"
    end
    
    return T, disturbance_list
end

# Build the transition matrix
const T, disturbance_list = build_transition_matrix()

##########################################################
# 3) Reward function
##########################################################
# +1 if in safe region, else 0
function reward(s::Tuple{Float64, Float64}, a::Float64)
    return is_safe(s[1], s[2]) ? 1.0 : 0.0
end

# Pre-compute rewards for all state-action pairs to avoid recomputation
function precompute_rewards()
    rewards = zeros(Float64, nstates, nactions)
    for s in 1:nstates
        for a in 1:nactions
            rewards[s, a] = reward(states_2d[s], actions[a])
        end
    end
    return rewards
end

const precomputed_rewards = precompute_rewards()

##########################################################
# 4) Solve the linear program (flow constraints)
##########################################################
function solve_case()
    println("*************************************************")
    println("Solving case with sigma = $sigma")
    println("*************************************************")
    
    # Set Gurobi license credentials right before creating the model
    ENV["GRB_WLSACCESSID"] = "52eb20bf-115c-42d3-931f-47561460611c"
    ENV["GRB_WLSSECRET"] = "587b8f93-6d53-43c9-af49-6b96ac589004"
    ENV["GRB_LICENSEID"] = "2611020"

    # Setup model
    model = Model(Gurobi.Optimizer)
    
    # Set solver parameters to improve performance
    set_optimizer_attribute(model, "Method", 2)  # Use barrier method (interior point)
    set_optimizer_attribute(model, "Threads", Threads.nthreads())  # Use all CPU threads
    set_optimizer_attribute(model, "OutputFlag", 1)  # Show solver output
    
    # Define z[s, a] and y[s, a]
    @variable(model, z[1:nstates, 1:nactions] >= 0)
    @variable(model, y[1:nstates, 1:nactions] >= 0)

    # c1: total measure of z is 1
    @constraint(model, c1, sum(z[s,a] for s in 1:nstates, a in 1:nactions) == 1)

    # c2: flow constraints - compute these more efficiently
    @constraint(model, c2[j=1:nstates],
        sum(z[s,a] * T[s,a,j] for s in 1:nstates, a in 1:nactions)
        == sum(z[j,a2] for a2 in 1:nactions)
    )

    # c3: alpha distribution
    # For demonstration, alpha = 1 / nstates (uniform)
    alpha_dist = 1.0 / nstates
    @constraint(model, c3[j=1:nstates],
        sum(z[j,a] for a in 1:nactions)
      + sum(y[j,a] for a in 1:nactions)
      - sum(y[s,a] * T[s,a,j] for s in 1:nstates, a in 1:nactions)
        == alpha_dist
    )

    # Objective: maximize sum_{s,a} z[s,a] * reward(s,a)
    # Use precomputed rewards for better performance
    @objective(model, Max, sum(z[s,a] * precomputed_rewards[s,a]
                              for s in 1:nstates, a in 1:nactions))

    # Solve the model
    optimize!(model)
    
    # Check the solution status
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
    Nx = length(theta)
    Ny = length(theta_dot)
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
    X = [theta[i] for j in 1:Ny, i in 1:Nx]  # Nx×Ny
    Y = [theta_dot[j] for j in 1:Ny, i in 1:Nx]

    # Retrieve duals for c2, c3
    dual_c2_vals = [dual(c2[j]) for j in 1:nstates]
    dual_c3_vals = [dual(c3[j]) for j in 1:nstates]

    # Reshape them into Nx×Ny
    G_map = zeros(Nx, Ny)
    H_map = zeros(Nx, Ny)
    for i in 1:Nx
        for j in 1:Ny
            sidx = (i - 1)*Ny + j
            G_map[i,j] = dual_c2_vals[sidx]
            H_map[i,j] = dual_c3_vals[sidx]
        end
    end

    return X, Y, Z_occup, G_map, H_map
end

##########################################################
# 5) Plot in 3D (z, g, h) side-by-side (No Colorbars)
#    Then save results to a new directory "results/"
##########################################################
function main_3D()
    # Solve the LP
    X, Y, Z_occup, G_map, H_map = solve_case()

    # "Safe set" boundary lines (for reference)
    safe_x = 0.3*[-1, 1, 1, -1, -1]
    safe_v = 0.6*[-1, -1, 1, 1, -1]
    safe_z = zeros(length(safe_x))
 
    # Create a figure with 3 subplots (no colorbars)
    fig = Figure(resolution=(1800, 600))

    # (1) Occupation measure z(s)
    ax1 = Axis3(fig[1,1],
        title  = "sum(value(z[s,a]) for a in 1:nactions",
        xlabel = "theta (θ)",
        ylabel = " theta_dot(θ̇)",
        zlabel = "Value"
    )
    # Makie uses row-major for surface, so pass Z_occup' (transpose)
    surface!(ax1, X, Y, Z_occup', colormap=:plasma)
    lines!(ax1, safe_x, safe_v, safe_z, color=:red, linewidth=3)

    # (2) Dual w.r.t. flow constraint: G_map
    ax2 = Axis3(fig[1,2],
        title  = "Dual Constraint 2",
        xlabel = "theta (θ)",
        ylabel = " theta_dot(θ̇)",
        zlabel = "Value"
    )
    surface!(ax2, X, Y, G_map', colormap=:viridis)
    lines!(ax2, safe_x, safe_v, safe_z, color=:red, linewidth=3)

    # (3) Dual w.r.t. alpha-dist constraint: H_map
    ax3 = Axis3(fig[1,3],
        title  = "Dual Constraint 3",
        xlabel = "theta (θ)",
        ylabel = " theta_dot(θ̇)",
        zlabel = "Value"
    )
    surface!(ax3, X, Y, H_map', colormap=:hot)
    lines!(ax3, safe_x, safe_v, safe_z, color=:red, linewidth=3)

    display(fig)
    println("\nAll subplots shown in one row (no colorbars).")

    # -------------------------------------------------------
    # Create a "results" directory, save the figure there.
    # -------------------------------------------------------
    results_dir = joinpath(@__DIR__, "results", string(sigma))
    mkpath(results_dir)  # Create the directory if it doesn't exist
    GLMakie.save(joinpath(results_dir, "plot.png"), fig)  # Save as a .png image
    println("Figure saved to $(joinpath(results_dir, "plot.png"))")

    return fig
end

##########################################################
# 6) Run
##########################################################
@time main_3D()