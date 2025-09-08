# Original code by Saber Omidi, 
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
# Set to 'false' to run the one-player benchmark (Fig 1 from the paper).
# Set to 'true' to run your two-player (control vs. disturbance) game.
const IS_TWO_PLAYER_GAME = true 

# ---------------------------
# System and grid parameters
# ---------------------------
const DT = 0.1
const DISCOUNT_RATE = 0.1 # This is λ
const GAMMA = exp(-DISCOUNT_RATE * DT)

# State grid
const X1_MIN, X1_MAX = -1.0, 5.0
const X2_MIN, X2_MAX = -5.0, 5.0
const NUM_POINTS_STATE_1 = 161
const NUM_POINTS_STATE_2 = 161

# Control grid
const U_MIN, U_MAX = -2.0, 2.0
const NUM_POINTS_ACTIONS = 41

# Value iteration parameters
const TOLERANCE = 1e-6 # Epsilon for convergence
const MAX_ITER = 1000

# --- Key Parameter from Paper ---
const TAU_BAR = 2.0

const L = sqrt(5)

# Safe Box K = [K1_MIN, K1_MAX] x [K2_MIN, K2_MAX]
const K1_MIN, K1_MAX = 0.0, 4.0
const K2_MIN, K2_MAX = -3.0, 3.0

# Disturbance parameters (only used if IS_TWO_PLAYER_GAME is true)
const SIGMA = 1.0
const MEAN = 0.0
const NUM_SAMPLES = 100
const d_min = -1.0
const d_max =  1.0 

# ---------------------------
# Signed distance function to Safe Set K
# ---------------------------
"""
Calculates the standard signed distance to a box.
Returns < 0 for points inside the box.
Returns > 0 for points outside the box.
"""
function signed_distance_to_box(px, py, box_x_min, box_x_max, box_y_min, box_y_max)
    # Distance from point to box exterior
    dx = max(box_x_min - px, 0.0, px - box_x_max)
    dy = max(box_y_min - py, 0.0, py - box_y_max)
    dist_outside = sqrt(dx^2 + dy^2)

    if dist_outside > 0.0
        return dist_outside
    else
        # Negative distance to the closest boundary for points inside
        dist_to_xmin = px - box_x_min
        dist_to_xmax = box_x_max - px
        dist_to_ymin = py - box_y_min
        dist_to_ymax = box_y_max - py
        return -min(dist_to_xmin, dist_to_xmax, dist_to_ymin, dist_to_ymax)
    end
end

# ---------------------------
# Double integrator dynamics
# ---------------------------
function di_dynamic(x, v, u, d)
    x_next = x + v * DT
    v_next = v + (u + d) * DT
    return [x_next, v_next]
end

# ---------------------------
# Surface function l(x) and h(x)
# ---------------------------
function compute_surface_functions(x1_grid, x2_grid, L)
    nx1, nx2 = length(x1_grid), length(x2_grid)
    l = zeros(nx1, nx2)

    for (i, xi) in enumerate(x1_grid)
        for (j, vj) in enumerate(x2_grid)
            s_K = signed_distance_to_box(xi, vj, K1_MIN, K1_MAX, K2_MIN, K2_MAX)
            l_val = -s_K
            l[i, j] = min(max(l_val, -L), L)
        end
    end

    h = l .- L
    return l, h
end

# ---------------------------
# Main script execution
# ---------------------------
function main()
    # --- Grid and Initialization ---
    x1_grid = collect(LinRange(X1_MIN, X1_MAX, NUM_POINTS_STATE_1))
    x2_grid = collect(LinRange(X2_MIN, X2_MAX, NUM_POINTS_STATE_2))
    
    state_2d = [(x, v) for x in x1_grid for v in x2_grid]
    nstates = length(state_2d)
    
    actions = collect(LinRange(U_MIN, U_MAX, NUM_POINTS_ACTIONS))
    println("Number of states: $nstates, Number of actions: $(length(actions))")

    # --- Disturbance setup ---
    D_list = [clamp(rand(Normal(MEAN, SIGMA)), d_min, d_max) for _ in 1:NUM_SAMPLES]
    
    if IS_TWO_PLAYER_GAME
        println("Running in TWO-PLAYER mode.")
    else
        println("Running in ONE-PLAYER (benchmark) mode.")
    end

    _, h_matrix = compute_surface_functions(x1_grid, x2_grid, L)
    
    # Using 'U' to align with paper's notation for the discounted value function
    U = vec(copy(h_matrix)) 
    h_vec = copy(U)
    U_next = similar(U)

    states_matrix = hcat([collect(s) for s in state_2d]...)
    tree = KDTree(states_matrix)

    println("\nStarting value iteration for λ = $DISCOUNT_RATE...")
    diff = Inf
    iteration = 0

    @time while diff > TOLERANCE && iteration < MAX_ITER
        iteration += 1
        max_diff = 0.0

        for state_index in 1:nstates
            x, v = state_2d[state_index]
            best_over_u = -Inf

            for u in actions
                if IS_TWO_PLAYER_GAME
                    # --- TWO-PLAYER LOGIC (Control vs. Disturbance) ---
                    worst_over_d = Inf
                    for d in D_list
                        next_state = di_dynamic(x, v, u, d)
                        idxs, _ = knn(tree, next_state, 1)
                        j = idxs[1]
                        val_at_neighbor = GAMMA * U[j]
                        worst_over_d = min(worst_over_d, val_at_neighbor)
                    end
                    best_over_u = max(best_over_u, worst_over_d)
                else
                    # --- ONE-PLAYER LOGIC (Control only) ---
                    next_state = di_dynamic(x, v, u, 0.0) # Disturbance is zero
                    idxs, _ = knn(tree, next_state, 1)
                    j = idxs[1]
                    val_at_neighbor = GAMMA * U[j]
                    best_over_u = max(best_over_u, val_at_neighbor)
                end
            end

            U_next[state_index] = min(h_vec[state_index], best_over_u)
            
            current_diff = abs(U_next[state_index] - U[state_index])
            if current_diff > max_diff
                max_diff = current_diff
            end
        end

        U .= U_next
        diff = max_diff

        if iteration % 5 == 0
            @printf "Iteration %4d, max diff = %.6f\n" iteration diff
        end
    end

    if diff <= TOLERANCE
        println("Value iteration converged in $iteration iterations.")
    else
        @warn "Value iteration did not converge within $MAX_ITER iterations."
    end
        
    # Calculate Z from the converged U
    Z = reshape(U .+ L, (NUM_POINTS_STATE_1, NUM_POINTS_STATE_2))

    # --- Save Results ---
    script_dir = @__DIR__
    results_dir = joinpath(script_dir, "results")
    if !isdir(results_dir)
        mkdir(results_dir)
    end
    
    output_path = joinpath(results_dir, "Value_function.txt")
    writedlm(output_path, Z, ',')
    println("Final value function is saved to: $output_path")

    # ----------------------------------------------------
    # --- NEW: Node Count Analysis and Save to File ---
    # ----------------------------------------------------
    println("Performing node count analysis...")

    # Calculate the Lower Bound (V function) from Z
    coef_lower_bound = exp(-DISCOUNT_RATE * TAU_BAR)
    Lower_bound_V = (Z .- L * (1 - coef_lower_bound)) ./ coef_lower_bound

    # Total number of nodes
    total_nodes = NUM_POINTS_STATE_1 * NUM_POINTS_STATE_2

    # Count for the Upper Bound (Z >= 0)
    nodes_in_upper_bound = sum(Z .>= 0)
    percent_upper = (nodes_in_upper_bound / total_nodes) * 100

    # Count for the Lower Bound (V >= 0)
    nodes_in_lower_bound = sum(Lower_bound_V .>= 0)
    percent_lower = (nodes_in_lower_bound / total_nodes) * 100

    # --- Write results to a text file ---
    analysis_path = joinpath(results_dir, "analysis_summary.txt")
    open(analysis_path, "w") do f
        write(f, "--- Node Count & Percentage Results ---\n\n")

        line1 = @sprintf("Upper Bound (Z>=0): %d / %d nodes (%.2f%%)\n",
                        nodes_in_upper_bound, total_nodes, percent_upper)
        write(f, line1)

        line2 = @sprintf("Lower Bound (V>=0): %d / %d nodes (%.2f%%)\n",
                        nodes_in_lower_bound, total_nodes, percent_lower)
        write(f, line2)
    end
    println("Analysis summary saved to: $analysis_path")
end

# Run the main function
main()
