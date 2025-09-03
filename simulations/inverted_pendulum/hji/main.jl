# Original code by Saber Omidi, 

using LinearAlgebra
using NearestNeighbors
using Printf
using DelimitedFiles
using Random
using Distributions

println("Packages imported.")

Random.seed!(123)
# ---------------------------
# System and grid parameters
# ---------------------------
const DT = 0.1
const DISCOUNT_RATE = 0.1 # This is λ
const GAMMA = exp(-DISCOUNT_RATE * DT)

# State grid
const X1_MIN, X1_MAX = -0.5, 0.5
const X2_MIN, X2_MAX = -1.0, 1.0
const NUM_POINTS_STATE_1 = 201
const NUM_POINTS_STATE_2 = 201

# Control grid
const U_MIN, U_MAX = -3.0, 3.0
const NUM_POINTS_ACTIONS = 21

# Value iteration parameters
const TOLERANCE = 1e-4 # Epsilon for convergence
const MAX_ITER = 1000

# --- Key Parameter from Paper ---
# τ̄ (tau_bar) is the assumed upper bound on the time-to-target.
# It is used to calculate the over-approximation boundary.
const TAU_BAR = 2.0

# Safe Box K = [K1_MIN, K1_MAX] x [K2_MIN, K2_MAX]
const K1_MIN, K1_MAX = -0.3, 0.3
const K2_MIN, K2_MAX = -0.6, 0.6

# Disturbance parameters for worst-case analysis
const SIGMA = 0.0
const MEAN = 0.0
const NUM_SAMPLES = 100
const d_min = -1.0
const d_max =  1.0

# ---------------------------
# Signed distance function to Safe Set K
# ---------------------------
"""
Calculates the signed distance to the safe set K (a box).
Returns > 0 outside the box, < 0 inside the box.
"""
function signed_distance_to_K(x, y)
    dx = max(K1_MIN - x, 0.0, x - K1_MAX)
    dy = max(K2_MIN - y, 0.0, y - K2_MAX)
    dist_outside = norm([dx, dy])

    if dist_outside > 0
        return dist_outside
    else
        dist_to_boundary = min(x - K1_MIN, K1_MAX - x, y - K2_MIN, K2_MAX - y)
        return -dist_to_boundary
    end
end


# ---------------------------
# Double integrator dynamics
# ---------------------------
function di_dynamic(x, v, u, d)
    #system parameters 
	m=2.0
	g=10.0
	l=1.0

    x_next = x + v * DT
    v_next = v + (g/l*sin(x)+1/m*l^2*u + d) * DT
    return [x_next, v_next]
end

# ---------------------------
# Surface function l(x) and h(x)
# ---------------------------
function compute_surface_functions(x1_grid, x2_grid)
    nx1, nx2 = length(x1_grid), length(x2_grid)
    l = zeros(nx1, nx2)
    for (i, xi) in enumerate(x1_grid)
        for (j, vj) in enumerate(x2_grid)
            l[i, j] = -signed_distance_to_K(xi, vj)
        end
    end

    global L
    L = maximum(abs.(l))
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
    
    state_2d = [(x, v) for x in x1_grid for v in x2_grid]
    nstates = length(state_2d)
    
    actions = collect(LinRange(U_MIN, U_MAX, NUM_POINTS_ACTIONS))
    println("Number of states: $nstates, Number of actions: $(length(actions))")

    # --- Disturbance setup ---
    D_list = [clamp(rand(Normal(MEAN, SIGMA)), d_min, d_max) for _ in 1:NUM_SAMPLES]
    println("Disturbance samples created with σ = $SIGMA.")

    min_val, max_val = extrema(D_list)
    println("Minimum: $min_val, Maximum: $max_val") # Output: Minimum: 1, Maximum: 9


    _, h_matrix, L = compute_surface_functions(x1_grid, x2_grid)

    V = vec(copy(h_matrix))
    h_vec = copy(V)
    V_next = similar(V)

    # --- KD-tree for nearest neighbor lookups ---
    states_matrix = hcat([collect(s) for s in state_2d]...)
    tree = KDTree(states_matrix)

    # --- Value Iteration ---
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
                worst_over_d = Inf
                for d in D_list
                    next_state = di_dynamic(x, v, u, d)
                    idxs, _ = knn(tree, next_state, 1) # Corrected function call
                    j = idxs[1]
                    val_at_neighbor = GAMMA * V[j]
                    worst_over_d = min(worst_over_d, val_at_neighbor)
                end
                best_over_u = max(best_over_u, worst_over_d)
            end

            V_next[state_index] = min(h_vec[state_index], best_over_u)
            
            current_diff = abs(V_next[state_index] - V[state_index])
            if current_diff > max_diff
                max_diff = current_diff
            end
        end

        V .= V_next
        diff = max_diff

        if iteration % 10 == 0
            @printf "Iteration %4d, max diff = %.6f\n" iteration diff
        end
    end

    if diff <= TOLERANCE
        println("Value iteration converged in $iteration iterations.")
    else
        @warn "Value iteration did not converge within $MAX_ITER iterations."
    end
        
    Z = reshape(V .+ L, (NUM_POINTS_STATE_1, NUM_POINTS_STATE_2))

    # --- Save Results for Server ---
    script_dir = @__DIR__
    results_dir = joinpath(script_dir, "results")
    if !isdir(results_dir)
        mkdir(results_dir)
        println("Created folder: $results_dir")
    end

    # Calculate levels for approximations based on the paper's theory
    under_approx_level = 0.0
    over_approx_level = L * (1 - exp(-DISCOUNT_RATE * TAU_BAR))
    
    @printf "\nUnder-approximation level (Z(x)=0): %.4f\n" under_approx_level
    @printf "Over-approximation level (using τ̄=%.1f): %.4f\n" TAU_BAR over_approx_level

    # Save everything to a single, self-contained text file for later plotting
    output_path = joinpath(results_dir, "value_function_and_levels_deter.txt")
    
    open(output_path, "w") do f
        # Write header with metadata for easy parsing later
        println(f, "# Results for Double Integrator MDR Simulation")
        println(f, "# Under-Approximation Level (Z=0)")
        writedlm(f, [under_approx_level], ',')
        println(f, "# Over-Approximation Level (Z=L(1-exp(-λτ̄)))")
        writedlm(f, [over_approx_level], ',')
        println(f, "# Z Value Function Matrix ($(NUM_POINTS_STATE_1)x$(NUM_POINTS_STATE_2))")
        
        # Write the Z matrix
        writedlm(f, Z, ',')
    end

    println("Final value function and approximation levels saved to: $output_path")

end

# Run the main function
main()