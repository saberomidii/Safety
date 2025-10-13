# Packages 
using Random, Distributions
using NearestNeighbors
using SparseArrays
using Mosek
using JuMP
import MathOptInterface as MOI
using MosekTools  

using Plots

using Dates
using DelimitedFiles
# Disturbance 
Random.seed!(10) # Setting the seed
μ=0.0
σ=1.0
nsamples = 100
d_list =rand(Normal(μ, σ), nsamples)
max_d = maximum(d_list)
min_d = minimum(d_list)
l_d = -1.0
up_d = 1.0
d_list_bounded = clamp.(d_list,l_d,up_d)
# Transition matrix 
num_points_state_1 = 81
num_points_state_2 = 81
num_points_state_3 = 21

x1 = collect(LinRange(-0.5, 0.5, num_points_state_1))
x2 = collect(LinRange(-1.0, 1.0, num_points_state_2))
x3 = collect(LinRange(0, 2π, num_points_state_3))

u  = collect(LinRange(-2.0, 2.0, 41))

const a_outer = 0.4
const b_outer = 0.9

const a_inner = 0.2
const b_inner = 0.70

const x_c = 0.0
const v_c = 0.0

function is_safe(x::Float64, v::Float64)
    # Center of both ellipses
    
    # Check if point is inside the larger ellipse
    inside_outer = ((x - x_c)^2) / a_outer^2 + ((v - v_c)^2) / b_outer^2 <= 1

    # Check if point is inside the smaller ellipse
    inside_inner = ((x - x_c)^2) / a_inner^2 + ((v - v_c)^2) / b_inner^2 <= 1

    # Safe set is the ring: inside the outer, but outside the inner
    return inside_outer && !inside_inner
end


# Continuous (noisy) dynamics for double-integrator
function di_dynamics(x1::Float64, x2::Float64, x3::Float64, u::Float64, d::Float64)    
    if !is_safe(x1, x2)
        return (x1, x2, x3)   # no movement if outside the safe region
    else
        dt=0.1
        V=0.1       #Constant Speed 
        x1_next = x1 + V*cos(x3)*dt
        x2_next = x2 + V*cos(x3)*dt
        x3_next = x3 + (u+d)*dt
        return (x1_next, x2_next, x3_next)
    end
end

# Cartesian product of (x, v) forms the entire state space
const states_3d = [(x, v, theta) for x in x1 for v in x2 for theta in x3]
const nstates = length(states_3d)

threshold_for_transit = 0.0
# Action discretization
nactions = length(u)

# Build a KD-tree for snapping continuous next-states to the nearest discrete state
states_matrix = hcat([collect(s) for s in states_3d]...)  # shape: 3 x nstates
tree = KDTree(states_matrix)

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
        for i in 1:nsamples
                (xn, vn, thetan) = di_dynamics(s[1], s[2],s[3], u[a], d_list_bounded[i])
                # For knn, pass a 2-element Vector, not a Tuple
                idxs, dists = knn(tree, [xn, vn, thetan], 1)
                # if first(idxs) != is || j > 1
		        T[is][ a, first(idxs)] += 1.0 / nsamples
        end
    end
 end
end

println("Negative Values checking") 
@assert all(t->all(t.>=0.0),T) "Negative Value!"

println("Sum is one checking") 
for action in 1:nactions
	@assert all(i -> sum(T[i][action,:])≈1.0, axes(T,1)) "Not All row sums are approximately 1!"
end

println("Removing small values")
for s in 1:nstates
	T[s]= dropzeros!(map(x->abs(x)<threshold_for_transit ? 0.0 : x, T[s]))
end

println("Done Building T.\n")
# AVR 
# Initialize reward vector for all states
# +1 if in safe region, else 0
function reward(s::Tuple{Float64, Float64, Float64})
    return is_safe(s[1], s[2]) ? 1.0 : 0.0
end


# Initialize reward vector for all states
r = zeros(nstates)
for s in 1:nstates
    r[s] = reward(states_3d[s])
end

# LP 
# Setup model
model = Model(Mosek.Optimizer)
# TOLERANCE SETTINGS FOR A LINEAR PROGRAM
# Setting this to 0 (MSK_OFF) disables the basis identification procedure.
set_optimizer_attribute(model, "MSK_IPAR_INTPNT_BASIS", 0)
# set_time_limit_sec(model, 1800.0)

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

# print_model_summary(model)
optimize!(model)

# Check solution status
stat = termination_status(model)
println("Solver status: $stat")

println("Solution found with status: $stat")
println("Objective value = ", objective_value(model))
            
# Retrieve the g and h values (optimal or best found)
g_opt = value.(g)
h_opt = value.(h)

# Check for NaN or Inf values and replace them
# g_opt = replace(g_opt, NaN => 0.0, Inf => 1.0, -Inf => 0.0)
# h_opt = replace(h_opt, NaN => 0.0, Inf => 1.0, -Inf => 0.0)
# Reshape for visualization
gain = reshape(g_opt, num_points_state_1, num_points_state_2,num_points_state_3)
bias = reshape(h_opt, num_points_state_1, num_points_state_2,num_points_state_3)
# Extract optimal policy (for each state, find action that maximizes the value)
optimal_policy = zeros(Int, nstates)
@time optimal_policy = [argmax(r[s] .+ (T[s] * h_opt)) for s in 1:nstates]
policy= reshape(optimal_policy,num_points_state_1,num_points_state_2,num_points_state_3)


sim_dir = joinpath(@__DIR__, "..", "Simulation")
# 2. Make sure the directory exists (create it if it doesn't)
mkpath(sim_dir)
# 3. Define the full path to the file
policy_path = joinpath(sim_dir, "optimal_policy_avr.csv")
writedlm(policy_path, policy, ';')
println("Optimal Policy saved to optimal_policy_avr.csv")

println("--- AVR is done")

# MDR 
function signed_distance_to_ellipse_ring(x1::Float64, x2::Float64)
    # Level set value for the outer ellipse.
    # This is < 0 inside, > 0 outside.
    dist_from_outer = ((x1 - x_c)^2) / a_outer^2 + ((x2 - v_c)^2) / b_outer^2 - 1.0

    dist_from_inner_hole = ((x1 - x_c)^2) / a_inner^2 + ((x2 - v_c)^2) / b_inner^2 - 1.0


    return max(dist_from_outer, -dist_from_inner_hole)
end

x_range = collect(LinRange(-2.0, 2.0, 101))
y_range = collect(LinRange(-2.0, 2.0, 101))

function VI_MDR(λ::Float64)
    max_iteration= 10000
    max_tolerance = 1e-9
    dt = 0.1
    γ= exp(-λ*dt)

    L = maximum(abs.([signed_distance_to_ellipse_ring(s[1], s[2]) for s in states_3d]))
    h = zeros(Float64, nstates)
    
    for s_idx in 1:nstates
        s = states_3d[s_idx]
        dist = signed_distance_to_ellipse_ring(s[1], s[2])
        l_val = -dist
        l_bounded = clamp(l_val, -L, L)
        h[s_idx] = l_bounded - L
    end

    println("Terminal cost vector h calculated.")
    @assert all(h .<= 1e-9) "Positive Value in h!"

    U = copy(h)
    U_prev = similar(U)
    mdr_policy = zeros(Int, nstates)
    iteration = 0
    diff = Inf

    while diff > max_tolerance && iteration < max_iteration
            iteration += 1
            copyto!(U_prev, U)

        # Iterate over the 1D state index
        for s_idx in 1:nstates
        best_val_over_actions = -Inf
        best_action_idx = 1

            for a_idx in 1:nactions
                # Use the original T format: a Vector of Sparse Matrices
                s_primes_indices, _ = findnz(T[s_idx][a_idx, :])

                # Worst-case logic operating on the 1D U vector
                min_val_over_next_states = if isempty(s_primes_indices)
                                                -Inf # No transition possible, worst outcome
                                           else
                # Find the minimum U value among all possible next states
                                                minimum(U[s_prime] for s_prime in s_primes_indices)
                                            end

                if min_val_over_next_states > best_val_over_actions
                    best_val_over_actions = min_val_over_next_states
                    best_action_idx = a_idx
                end
            end

            U[s_idx] = min(h[s_idx], γ * best_val_over_actions)
            mdr_policy[s_idx] = best_action_idx
        end

        diff = maximum(abs.(U .- U_prev))
        if iteration % 1 == 0
             println("Iteration: ", iteration, "  Diff:", diff)
        end
    end

    println("Value Iteration loop is over after $iteration iterations.")
    # --- PART 4: Final Calculation and Plotting ---
    Z = U .+ L
    # Reshape the final 1D Z vector into a 2D map for plotting
    Z_map = reshape(Z, num_points_state_1, num_points_state_2, num_points_state_3) # Transpose to match grid orientation

    optimal_policy = zeros(Int, nstates)
    for s_idx in 1:nstates
        best_action_value = -Inf
        best_action_idx = 1
        for a_idx in 1:nactions
            s_primes_indices, _ = findnz(T[s_idx][a_idx, :])
            current_action_value = if isempty(s_primes_indices)
                -Inf
            else
                # This line will now work correctly because U_vec_00 is a vector
                minimum(U[s_prime] for s_prime in s_primes_indices)
            end
            if current_action_value > best_action_value
                best_action_value = current_action_value
                best_action_idx = a_idx
            end
        end
        optimal_policy[s_idx] = best_action_idx
    end

    policy_map = reshape(optimal_policy, num_points_state_1, num_points_state_2, num_points_state_3)


    return Z_map, Z, U, policy_map
end

λ= [0.0,0.1,0.2,0.3]
dt =0.1
MDR_MAX_ITERATION =10000
MDR_MAX_TOLERANCE = 1e-9
# --- Call the solver for each case ---
Z_map_00,Z_0,U_0,policy_0 = VI_MDR(λ[1])

Z_map_01,Z_1,U_1,policy_1 = VI_MDR(λ[2])

Z_map_02,Z_2,U_2,policy_2 = VI_MDR(λ[3])

Z_map_03,Z_3,U_3,policy_3 = VI_MDR(λ[4])

# --- Generate the multi-layered plot ---

# Create the initial plot and specify the legend's position
#p = contour(x1, x2, gain, levels=[1.0], color=:black, linewidth=3, label="AVR Safe Set (g=1)", legend=:outertopright)

# Add the MDR results using contour!
#contour!(p, x1, x2, Z_map_00, levels=[0.0], color=:red, linestyle=:dash, linewidth=2, label="MDR (λ=0.0)")
#contour!(p, x1, x2, Z_map_01, levels=[0.0], color=:purple, linestyle=:dashdot, linewidth=2, label="MDR (λ=0.01)")
#contour!(p, x1, x2, Z_map_02, levels=[0.0], color=:blue, linestyle=:dot, linewidth=2, label="MDR (λ=0.1)")
#contour!(p, x1, x2, Z_map_03, levels=[0.0], color=:green, linewidth=2, label="MDR (λ=0.2)")

# # Add labels and title
# xlabel!("Position (x1)")
# ylabel!("Velocity (x2)")
# title!("AVR Safe Set vs. MDR Reachable Set Boundaries")

# # Display the plot
# p

# ... after the VI_MDR calls ...

# --- Calculate the percentage of states where Z is approximately 0 ---
# Define a small tolerance for checking if a value is close to zero
zero_tolerance = 1e-9
z0_percentage_00 = 100 * count(z -> isapprox(z, 0.0, atol=zero_tolerance), Z_0) / nstates
z0_percentage_01 = 100 * count(z -> isapprox(z, 0.0, atol=zero_tolerance), Z_1) / nstates
z0_percentage_02 = 100 * count(z -> isapprox(z, 0.0, atol=zero_tolerance), Z_2) / nstates
z0_percentage_03 = 100 * count(z -> isapprox(z, 0.0, atol=zero_tolerance), Z_3) / nstates

# ==========================================================
# SAVE RESULTS
# ==========================================================
println("\n--- Saving results to file ---")

# --- 1. Create a Timestamped Directory for the Results ---
timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
results_dir = joinpath(@__DIR__, "results", timestamp)
mkpath(results_dir)
println("Results will be saved in: ", results_dir)


# --- 2. Write the Updated Summary File ---
summary_filepath = joinpath(results_dir, "summary.txt")
open(summary_filepath, "w") do f
    println(f, "Safety Analysis Results")
    println(f, "========================================")
    println(f, "Run Timestamp: ", timestamp)
    println(f, "System Type: Double Integrator")

    println(f, "\n--- System Details ---")
    println(f, "dt: ", dt)
    println(f, "safe_region: (x_lims = ($a_outer, $b_outer), v_lims = ($a_inner, $b_inner))")

    println(f, "\n--- Discretization ---")
    println(f, "Number of States: ", nstates)
    println(f, "Number of Actions: ", nactions)
    println(f, "State Grid (x1): ", num_points_state_1, " points from [", minimum(x1), ", ", maximum(x1), "]")
    println(f, "State Grid (x2): ", num_points_state_2, " points from [", minimum(x2), ", ", maximum(x2), "]")
    println(f, "State Grid (x3): ", num_points_state_3, " points from [", minimum(x3), ", ", maximum(x3), "]")

    println(f, "\n--- Disturbance Properties ---")
    println(f, "Distribution Type: Normal")
    println(f, "Distribution Parameters: (mean = $μ, sigma = $σ)")
    println(f, "Clamping Bounds: [", l_d, ", ", up_d, "]")
    println(f, "Samples Generated: ", nsamples)

    println(f, "\n--- Transition Matrix ---")
    println(f, "Monte Carlo Samples per (s,a): ", nstates)

    println(f, "\n--- Average Reward (AVR) Results ---")
    println(f, "Objective (min average reward): ", objective_value(model))
    g1_count = count(x -> isapprox(x, 1.0; atol=1e-3), g_opt)
    println(f, "Percentage of states with g(x) ≈ 1: ", round(100 * g1_count / nstates, digits=2), "%")

    println(f, "\n--- Minimum Discounted Reward (MDR) Results ---")
    println(f, "Note: Using worst-case (minimum) value iteration.")
    # --- NEW LINES START HERE ---
    println(f, "Percentage of states with Z(x) ≈ 0 (λ=$(λ[1])): ", "%")
    println(f, "Percentage of states with Z(x) ≈ 0 (λ=$(λ[2])): ",  "%")
    println(f, "Percentage of states with Z(x) ≈ 0 (λ=$(λ[3])): ", "%")
    println(f, "Percentage of states with Z(x) ≈ 0 (λ=$(λ[4])): ", "%")
    # --- NEW LINES END HERE ---

    println(f, "\nMDR Parameters:")
    println(f, "MAX_ITER: ", MDR_MAX_ITERATION)
    println(f, "TOLERANCE: ", MDR_MAX_TOLERANCE)
    println(f, "λ values tested: ", λ)
end
println("Summary file saved to: ", summary_filepath)

# --- 3. Save Data Matrices to CSV Files ---
writedlm(joinpath(results_dir, "AVR_gain_map.csv"), gain, ',')
writedlm(joinpath(results_dir, "MDR_Z_map_lambda_0.0.csv"), Z_map_00, ',')
writedlm(joinpath(results_dir, "MDR_Z_map_lambda_0.01.csv"), Z_map_01, ',')
writedlm(joinpath(results_dir, "MDR_Z_map_lambda_0.1.csv"), Z_map_02, ',')
writedlm(joinpath(results_dir, "MDR_Z_map_lambda_0.2.csv"), Z_map_03, ',')

writedlm(joinpath(sim_dir, "MDR_policy_map_lambda_0.0.csv"), policy_0, ',')
writedlm(joinpath(sim_dir, "MDR_policy_map_lambda_0.01.csv"), policy_1, ',')
writedlm(joinpath(sim_dir, "MDR_policy_map_lambda_0.1.csv"), policy_2, ',')
writedlm(joinpath(sim_dir, "MDR_policy_map_lambda_0.2.csv"), policy_3, ',')

println("All result matrices saved to CSV files.")
