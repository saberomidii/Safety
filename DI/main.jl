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
l_d =-1.0
up_d =1.0
d_list_bounded = clamp.(d_list,l_d,up_d)
# Transition matrix 
num_points_state_1 = 11
num_points_state_2 = 11
x1 = collect(LinRange(-1.0, 5.0, num_points_state_1))
x2 = collect(LinRange(-5.0, 5.0, num_points_state_2))
u  = collect(LinRange(-2.0, 2.0, 11))

const c_min_1 = 0.0
const c_max_1 = 4.0

const c_min_2 = -3.0
const c_max_2 =  3.0


function is_safe(x1::Float64, x2::Float64)
    # "Safe" region: x in [0,4], v in [-3,3]
    return (c_min_1 <= x1 <= c_max_1) && (c_min_2 <= x2 <= c_max_2)
end


# Continuous (noisy) dynamics for double-integrator
function di_dynamics(x1::Float64, x2::Float64, u::Float64, d::Float64)    
    if !is_safe(x1, x2)
        return (x1, x2)   # no movement if outside the safe region
    else
    dt=0.1
    x1_next = x1 + x2*dt
	x2_next = x2 + (u + d)*dt
        return (x1_next, x2_next)
    end
end


# Cartesian product of (x, v) forms the entire state space
states_2d = [(x, v) for x in x1 for v in x2]
nstates = length(states_2d)
threshold_for_transit = 0.0
# Action discretization
nactions = length(u)

# Build a KD-tree for snapping continuous next-states to the nearest discrete state
states_matrix = hcat([collect(s) for s in states_2d]...)  # shape: 3 x nstates
tree = KDTree(states_matrix)

# Initialize T array: T[s, a, s_next]
T = Vector{SparseMatrixCSC{Float64,Int64}}(undef,nstates)

for index_state in 1:nstates
	T[index_state] = spzeros(nactions,nstates)
end

# (Compute Transition Matrix)
@time begin
 for is in 1:nstates
    s = states_2d[is]   # e.g., s = (x, y, theta)
    for a in 1:nactions
        for i in 1:nsamples
                (xn, vn) = di_dynamics(s[1], s[2], u[a], d_list_bounded[i])
                # For knn, pass a 2-element Vector, not a Tuple
                idxs, dists = knn(tree, [xn, vn], 1)
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
function reward(s::Tuple{Float64, Float64})
    return is_safe(s[1], s[2]) ? 1.0 : 0.0
end


# Initialize reward vector for all states
r = zeros(nstates)
for s in 1:nstates
    r[s] = reward(states_2d[s])
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
println("\nStarting AVR LP optimization...")
# --- Start Timer for LP Solve ---
avr_lp_start_time = time()
optimize!(model)
avr_lp_solve_time = time() - avr_lp_start_time
# --------------------------------
println("AVR LP Solve Time: $(avr_lp_solve_time) seconds")

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
gain = reshape(g_opt, num_points_state_1, num_points_state_2)
bias = reshape(h_opt, num_points_state_1, num_points_state_2)
# Extract optimal policy (for each state, find action that maximizes the value)
optimal_policy = zeros(Int, nstates)
@time optimal_policy = [argmax(r[s] .+ (T[s] * h_opt)) for s in 1:nstates]
policy= reshape(optimal_policy,num_points_state_1,num_points_state_2)

println("--- AVR is done")

# MDR 

function signed_distance_to_box(x1, x2, c_min_1, c_max_1, c_min_2, c_max_2)
    dx = max(c_min_1 - x1, 0, x1 - c_max_1)
    dy = max(c_min_2 - x2, 0, x2 - c_max_2)
    dist_outside = sqrt(dx^2 + dy^2)
    return dist_outside > 0 ? dist_outside : -min(x1 - c_min_1, c_max_1 - x1, x2 - c_min_2, c_max_2 - x2)
end

function VI_MDR(λ::Float64)
    max_iteration= 10000
    max_tolerance = 1e-9
    dt = 0.1
    γ= exp(-λ*dt)

    L = maximum(abs.([signed_distance_to_box(s[1], s[2], c_min_1, c_max_1, c_min_2, c_max_2) for s in states_2d]))
    h = zeros(Float64, nstates)
    
    for s_idx in 1:nstates
        s = states_2d[s_idx]
        dist = signed_distance_to_box(s[1], s[2], c_min_1, c_max_1, c_min_2, c_max_2)
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

    # --- Start Timer for the Value Iteration Loop ---
    vi_start_time = time()

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

    vi_solve_time = time() - vi_start_time
    println("Value Iteration loop is over after $iteration iterations.")
    # --- PART 4: Final Calculation and Plotting ---
    Z = U .+ L
    # Reshape the final 1D Z vector into a 2D map for plotting
    Z_map = reshape(Z, num_points_state_1, num_points_state_2) # Transpose to match grid orientation



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

    policy_map = reshape(optimal_policy, num_points_state_1, num_points_state_2)'


    return Z_map, Z, U, policy_map, vi_solve_time
end

λ= [0.0,0.01,0.015,0.02]
dt =0.1
MDR_MAX_ITERATION =10000
MDR_MAX_TOLERANCE = 1e-9

# --- Call the solver for each case and capture the time ---
Z_map_00,Z_0,U_0,policy_0, time_00 = VI_MDR(λ[1])
Z_map_01,Z_1,U_1,policy_1, time_01 = VI_MDR(λ[2])
Z_map_02,Z_2,U_2,policy_2, time_02 = VI_MDR(λ[3])
Z_map_03,Z_3,U_3,policy_3, time_03 = VI_MDR(λ[4])



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

# Make sure you have these at the top of your script
using Dates
using DelimitedFiles

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
    println(f, "safe_region: (x_lims = ($c_min_1, $c_max_1), v_lims = ($c_min_2, $c_max_2))")

    println(f, "\n--- Discretization ---")
    println(f, "Number of States: ", nstates)
    println(f, "Number of Actions: ", nactions)
    println(f, "State Grid (x1): ", num_points_state_1, " points from [", minimum(x1), ", ", maximum(x1), "]")
    println(f, "State Grid (x2): ", num_points_state_2, " points from [", minimum(x2), ", ", maximum(x2), "]")

    println(f, "\n--- Disturbance Properties ---")
    println(f, "Distribution Type: Normal")
    println(f, "Distribution Parameters: (mean = $μ, sigma = $σ)")
    println(f, "Clamping Bounds: [", l_d, ", ", up_d, "]")
    println(f, "Samples Generated: ", nsamples)

    println(f, "\n--- Transition Matrix ---")
    println(f, "Monte Carlo Samples per (s,a): ", nstates)

    println(f, "\n--- Average Reward (AVR) Results ---")
    println(f, "LP Solve Time: $(avr_lp_solve_time) seconds")
    println(f, "Objective (min average reward): ", objective_value(model))
    g1_count = count(x -> isapprox(x, 1.0; atol=1e-3), g_opt)
    println(f, "Percentage of states with g(x) ≈ 1: ", round(100 * g1_count / nstates, digits=2), "%")

    println(f, "\n--- Minimum Discounted Reward (MDR) Results ---")
    println(f, "Note: Using worst-case (minimum) value iteration.")

    # --- ADD VI SOLVE TIMES HERE ---
    println(f, "MDR VI Solve Time (λ=$(λ[1])): $(time_00) seconds")
    println(f, "MDR VI Solve Time (λ=$(λ[2])): $(time_01) seconds")
    println(f, "MDR VI Solve Time (λ=$(λ[3])): $(time_02) seconds")
    println(f, "MDR VI Solve Time (λ=$(λ[4])): $(time_03) seconds")
    # -------------------------------
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
writedlm(joinpath(results_dir, "MDR_policy_map_lambda_0.0.csv"), policy_0, ',')
writedlm(joinpath(results_dir, "MDR_policy_map_lambda_0.01.csv"), policy_1, ',')
writedlm(joinpath(results_dir, "MDR_policy_map_lambda_0.1.csv"), policy_2, ',')
writedlm(joinpath(results_dir, "MDR_policy_map_lambda_0.2.csv"), policy_3, ',')

println("All result matrices saved to CSV files.")
