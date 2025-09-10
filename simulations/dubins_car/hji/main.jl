# Updated code based on original by Saber Omidi
# This script adapts the value iteration safety game from an inverted pendulum
# to a Dubins car model.

using LinearAlgebra
using NearestNeighbors
using Printf
using DelimitedFiles
using Random
using Distributions

println("Packages imported.")

Random.seed!(10)

# ---------------------------
# Game Setup Switch
# ---------------------------
# Set to 'false' to run the one-player benchmark.
# Set to 'true' to run a two-player (control vs. disturbance) game.
const IS_TWO_PLAYER_GAME = true

# ---------------------------
# System and grid parameters for Dubins Car
# ---------------------------
const DT = 0.1
const DISCOUNT_RATE = 0.1 # This is λ
const GAMMA = exp(-DISCOUNT_RATE * DT)

# State grid (x, y, theta)
const X1_MIN, X1_MAX = -2.0, 2.0
const X2_MIN, X2_MAX = -2.0, 2.0
const X3_MIN, X3_MAX = 0.0, 2π
const NUM_POINTS_STATE_1 = 101
const NUM_POINTS_STATE_2 = 101
const NUM_POINTS_STATE_3 = 18

# Control grid (steering rate)
const U_MIN, U_MAX = -1.0, 1.0 # Turn rate in rad/s
const NUM_POINTS_ACTIONS = 21

# Constant velocity
const V = 1.0 # m/s

# Value iteration parameters
const TOLERANCE = 1e-6
const MAX_ITER = 1000

# --- Key Parameter from Paper ---
const TAU_BAR = 2.0

# Safe Box K for Dubins Car (in x, y)
const K1_MIN, K1_MAX = -1, 1
const K2_MIN, K2_MAX = -2, 2
# The heading (theta) is not restricted for the safe set.

# Disturbance parameters (affecting steering)
const SIGMA_D = 1.0
const MEAN_D = 0.0
const NUM_SAMPLES = 100
const d_min = -1
const d_max =  1

# ---------------------------
# Signed distance function to Safe Set K (2D box)
# ---------------------------
function signed_distance_to_box_2d(px, py, box_x_min, box_x_max, box_y_min, box_y_max)
    dx = max(box_x_min - px, 0.0, px - box_x_max)
    dy = max(box_y_min - py, 0.0, py - box_y_max)
    dist_outside = sqrt(dx^2 + dy^2)

    if dist_outside > 0.0
        return dist_outside
    else
        dist_to_xmin = px - box_x_min
        dist_to_xmax = box_x_max - px
        dist_to_ymin = py - box_y_min
        dist_to_ymax = box_y_max - py
        return -min(dist_to_xmin, dist_to_xmax, dist_to_ymin, dist_to_ymax)
    end
end

# ---------------------------
# Dubins Car Dynamics
# ---------------------------
function dubins_dynamic(x, y, theta, u, d)
    # The control 'u' is the steering rate, 'd' is the disturbance.
    # The car's velocity 'V' is constant.
    
    x_next = x + V * cos(theta) * DT
    y_next = y + V * sin(theta) * DT
    theta_next = theta + (u + d) * DT
    
    return [x_next, y_next, theta_next]
end

# ---------------------------
# Surface function l(x) and h(x)
# ---------------------------
function compute_surface_functions(x1_grid, x2_grid, x3_grid)
    nx1, nx2, nx3 = length(x1_grid), length(x2_grid), length(x3_grid)
    
    # First, calculate the unclipped cost to find the maximum magnitude L
    l_unclipped = zeros(nx1, nx2, nx3)
    for (i, xi) in enumerate(x1_grid)
        for (j, yj) in enumerate(x2_grid)
            for (k, _) in enumerate(x3_grid)
                s_K = signed_distance_to_box_2d(xi, yj, K1_MIN, K1_MAX, K2_MIN, K2_MAX)
                l_unclipped[i, j, k] = -s_K
            end
        end
    end
    L = maximum(abs.(l_unclipped))
    println("Dynamically calculated L = $L")

    # Now, compute the final, clipped l(x) as specified in the paper
    l = zeros(nx1, nx2, nx3)
    for i in 1:nx1
        for j in 1:nx2
            for k in 1:nx3
                l_val = l_unclipped[i, j, k]
                l[i, j, k] = min(max(l_val, -L), L)
            end
        end
    end
    
    h = l .- L
    return l, h, L
end

# ---------------------------
# Main script execution
# ---------------------------
function main()
    # --- Grid and Initialization ---
    x1_grid = collect(LinRange(X1_MIN, X1_MAX, NUM_POINTS_STATE_1))
    x2_grid = collect(LinRange(X2_MIN, X2_MAX, NUM_POINTS_STATE_2))
    x3_grid = collect(LinRange(X3_MIN, X3_MAX, NUM_POINTS_STATE_3))
    
    # Create a list of all states
    state_3d = [(x, y, theta) for x in x1_grid for y in x2_grid for theta in x3_grid]
    nstates = length(state_3d)
    
    actions = collect(LinRange(U_MIN, U_MAX, NUM_POINTS_ACTIONS))
    println("Number of states: $nstates, Number of actions: $(length(actions))")

    D_list = [clamp(rand(Normal(MEAN_D, SIGMA_D)), d_min, d_max) for _ in 1:NUM_SAMPLES]
    if IS_TWO_PLAYER_GAME println("Running in TWO-PLAYER mode.") else println("Running in ONE-PLAYER (benchmark) mode.") end

    _, h_matrix, L = compute_surface_functions(x1_grid, x2_grid, x3_grid)
    
    # Vectorize h_matrix and initialize U
    U = vec(copy(h_matrix))
    h_vec = copy(U)
    U_next = similar(U)

    # Convert state list to matrix for KDTree
    states_matrix = hcat([collect(s) for s in state_3d]...)
    tree = KDTree(states_matrix)

    println("\nStarting value iteration for λ = $DISCOUNT_RATE...")
    diff = Inf
    iteration = 0

    @time while diff > TOLERANCE && iteration < MAX_ITER
        iteration += 1
        max_diff = 0.0
        for state_index in 1:nstates
            x, y, theta = state_3d[state_index]
            best_over_u = -Inf
            for u in actions
                if IS_TWO_PLAYER_GAME
                    worst_over_d = Inf
                    for d in D_list
                        next_state = dubins_dynamic(x, y, theta, u, d)
                        idxs, _ = knn(tree, next_state, 1)
                        j = idxs[1]
                        val_at_neighbor = GAMMA * U[j]
                        worst_over_d = min(worst_over_d, val_at_neighbor)
                    end
                    best_over_u = max(best_over_u, worst_over_d)
                else
                    next_state = dubins_dynamic(x, y, theta, u, 0.0)
                    idxs, _ = knn(tree, next_state, 1)
                    j = idxs[1]
                    val_at_neighbor = GAMMA * U[j]
                    best_over_u = max(best_over_u, val_at_neighbor)
                end
            end
            U_next[state_index] = min(h_vec[state_index], best_over_u)
            current_diff = abs(U_next[state_index] - U[state_index])
            if current_diff > max_diff max_diff = current_diff end
        end
        U .= U_next
        diff = max_diff
        if iteration % 10 == 0 @printf "Iteration %4d, max diff = %.6f\n" iteration diff end
    end

    if diff <= TOLERANCE println("Value iteration converged in $iteration iterations.") else @warn "Did not converge." end
        
    # Reshape U back into Z matrix.
    Z = reshape(U .+ L, (NUM_POINTS_STATE_1, NUM_POINTS_STATE_2, NUM_POINTS_STATE_3))

    # --- Save Results ---
    script_dir = @__DIR__
    results_dir = joinpath(script_dir, "results")
    if !isdir(results_dir) mkdir(results_dir) end
    
    output_path = joinpath(results_dir, "Value_function_dubins.txt")
    # Save the flattened Z matrix and a header for reshaping
    open(output_path, "w") do f
        write(f, "# Dimensions: $(NUM_POINTS_STATE_1),$(NUM_POINTS_STATE_2),$(NUM_POINTS_STATE_3)\n")
        writedlm(f, vec(Z), ',')
    end
    println("Final value function saved to: $output_path")

    # --- Node Count Analysis and Save to File ---
    println("Performing node count analysis...")
    coef_lower_bound = exp(-DISCOUNT_RATE * TAU_BAR)
    Lower_bound_V = (vec(Z) .- L * (1 - coef_lower_bound)) ./ coef_lower_bound
    total_nodes = NUM_POINTS_STATE_1 * NUM_POINTS_STATE_2 * NUM_POINTS_STATE_3
    nodes_in_upper_bound = sum(vec(Z) .>= 0)
    percent_upper = (nodes_in_upper_bound / total_nodes) * 100
    nodes_in_lower_bound = sum(Lower_bound_V .>= 0)
    percent_lower = (nodes_in_lower_bound / total_nodes) * 100

    analysis_path = joinpath(results_dir, "dubins_analysis.txt")
    open(analysis_path, "w") do f
        write(f, "--- Dubins Car Node Count & Percentage Results ---\n\n")
        line1 = @sprintf("Upper Bound (Z>=0): %d / %d nodes (%.2f%%)\n",
                         nodes_in_upper_bound, total_nodes, percent_upper)
        write(f, line1)
        line2 = @sprintf("Lower Bound (V>=0): %d / %d nodes (%.2f%%)\n",
                         nodes_in_lower_bound, total_nodes, percent_lower)
        write(f, line2)
    end
    println("Analysis summary saved to: $analysis_path")
end

main()


