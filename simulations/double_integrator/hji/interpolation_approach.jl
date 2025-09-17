# --- Required Packages ---
using Plots
using Printf
using LinearAlgebra
using DelimitedFiles # Add this line for saving text files

# ---------------------------
# GLOBAL CONSTANTS
# ---------------------------
const X1_MIN, X1_MAX = -1.0, 5.0
const X2_MIN, X2_MAX = -5.0, 5.0
const NUM_POINTS_STATE_1 = 161
const NUM_POINTS_STATE_2 = 161

const U_MAX = 2.0
const NUM_POINTS_ACTIONS = 41

const K1_MIN, K1_MAX = 0.0, 4.0
const K2_MIN, K2_MAX = -3.0, 3.0

const LAMBDA = 0.0
const TAU_BAR = 2.0
const DT = 0.1
const GAMMA = exp(-LAMBDA * DT)

const TOLERANCE = 1e-6
const MAX_ITER = 1000
const L = sqrt((5.0 - 4.0)^2 + (5.0 - 3.0)^2)

# ---------------------------
# HELPER FUNCTIONS
# ---------------------------
function signed_distance_to_box(px, py)
    dx = max(K1_MIN - px, 0.0, px - K1_MAX)
    dy = max(K2_MIN - py, 0.0, py - K2_MAX)
    dist_outside = sqrt(dx^2 + dy^2)
    if dist_outside > 0.0
        return dist_outside
    else
        dist_to_xmin = px - K1_MIN; dist_to_xmax = K1_MAX - px
        dist_to_ymin = py - K2_MIN; dist_to_ymax = K2_MAX - py
        return -min(dist_to_xmin, dist_to_xmax, dist_to_ymin, dist_to_ymax)
    end
end

function dynamics(x1, x2, u)
    x1_next = x1 + x2 * DT
    x2_next = x2 + u * DT
    return (x1_next, x2_next)
end

function interpolate_bilinear(x1p, x2p, x1_grid, x2_grid, U_matrix)
    x1p = clamp(x1p, X1_MIN, X1_MAX)
    x2p = clamp(x2p, X2_MIN, X2_MAX)
    i = searchsortedfirst(x1_grid, x1p) - 1
    j = searchsortedfirst(x2_grid, x2p) - 1
    i = max(1, min(i, NUM_POINTS_STATE_1 - 1))
    j = max(1, min(j, NUM_POINTS_STATE_2 - 1))
    x1_i, x1_i1 = x1_grid[i], x1_grid[i+1]
    x2_j, x2_j1 = x2_grid[j], x2_grid[j+1]
    U_ij  = U_matrix[i, j]
    U_i1j = U_matrix[i+1, j]
    U_ij1 = U_matrix[i, j+1]
    U_i1j1= U_matrix[i+1, j+1]
    tx = (x1p - x1_i) / (x1_i1 - x1_i)
    ty = (x2p - x2_j) / (x2_j1 - x2_j)
    U_bottom = (1 - tx) * U_ij + tx * U_i1j
    U_top    = (1 - tx) * U_ij1 + tx * U_i1j1
    return (1 - ty) * U_bottom + ty * U_top
end

# ---------------------------
# MAIN SIMULATION FUNCTION
# ---------------------------
function solve_and_plot()
    println("--- Starting MDR Value Iteration for Double Integrator ---")

    # --- Grid and Initialization ---
    x1_grid = collect(LinRange(X1_MIN, X1_MAX, NUM_POINTS_STATE_1))
    x2_grid = collect(LinRange(X2_MIN, X2_MAX, NUM_POINTS_STATE_2))
    actions = collect(LinRange(-U_MAX, U_MAX, NUM_POINTS_ACTIONS))
    nstates = NUM_POINTS_STATE_1 * NUM_POINTS_STATE_2
    println("Grid size: $NUM_POINTS_STATE_1 x $NUM_POINTS_STATE_2 = $nstates states.")
    h_matrix = [-signed_distance_to_box(x, v) - L for x in x1_grid, v in x2_grid]
    h_vec = vec(h_matrix)
    U = copy(h_vec)
    U_next = similar(U)

    # --- Value Iteration Loop ---
    println("\nStarting value iteration for λ = $LAMBDA...")
    iteration = 0
    diff = Inf
    @time while diff > TOLERANCE && iteration < MAX_ITER
        iteration += 1
        max_diff = 0.0
        U_matrix = reshape(U, (NUM_POINTS_STATE_1, NUM_POINTS_STATE_2))
        for i in 1:nstates
            ix1 = (i-1) % NUM_POINTS_STATE_1 + 1
            ix2 = (i-1) ÷ NUM_POINTS_STATE_1 + 1
            x1, x2 = x1_grid[ix1], x2_grid[ix2]
            best_over_u = -Inf
            for u in actions
                x1_next, x2_next = dynamics(x1, x2, u)
                val_at_neighbor = GAMMA * interpolate_bilinear(x1_next, x2_next, x1_grid, x2_grid, U_matrix)
                best_over_u = max(best_over_u, val_at_neighbor)
            end
            U_next[i] = min(h_vec[i], best_over_u)
            max_diff = max(max_diff, abs(U_next[i] - U[i]))
        end
        U .= U_next
        diff = max_diff
        if iteration % 10 == 0
            @printf "Iteration %4d, max diff = %.8f\n" iteration diff
        end
    end

    if diff <= TOLERANCE
        println("Value iteration converged in $iteration iterations.")
    else
        @warn "Value iteration did not converge within $MAX_ITER iterations."
    end
    
    # --- Post-Processing ---
    Z = reshape(U, (NUM_POINTS_STATE_1, NUM_POINTS_STATE_2)) .+ L

    # #######################################################################
    # ## NEW SECTION: Save Results to Text File
    # #######################################################################
    results_dir = joinpath(@__DIR__, "interpolation_results")
    mkpath(results_dir) # Creates the directory if it doesn't exist
    output_path = joinpath(results_dir, "Value_function_lambda_0.0.txt")
    println("\nSaving Z matrix to: $output_path")
    writedlm(output_path, Z, ',')
    println("--- Data saved successfully. ---")
    # #######################################################################

    # --- Plotting ---
    level_over = 0.0
    level_under = L * (1 - exp(-LAMBDA * TAU_BAR))
    @printf "\nPlotting Over-approx (Z(x)>%.2f) and Under-approx (Z(x)>%.4f)\n" level_over level_under

    safe_box_shape = Shape([K1_MIN, K1_MAX, K1_MAX, K1_MIN], [K2_MIN, K2_MIN, K2_MAX, K2_MAX])
    p = plot(safe_box_shape, opacity=0.3, color=:red, label="Safe Box K",
        xlabel="State x₁", ylabel="State x₂",
        title="MDR Safe Set Approximations (λ=$LAMBDA)",
        aspect_ratio=:equal,
        xlims=(X1_MIN, X1_MAX), ylims=(X2_MIN, X2_MAX),
        legend=:bottomleft,
        framestyle=:box
    )
    contour!(p, x1_grid, x2_grid, Z', levels=[level_over],
        linewidth=3, color=:blue, linestyle=:solid, label="Z(x) Over-Approximation")
    contour!(p, x1_grid, x2_grid, Z', levels=[level_under],
        linewidth=3, color=:darkblue, linestyle=:dash, label="Z(x) Under-Approximation")

    display(p)
    println("--- Plot generated successfully. ---")
end

# ---------------------------
# RUN THE SCRIPT
# ---------------------------
solve_and_plot()


# -----------------------------------------------------------------------------
# Complete Julia Script to Reproduce Figure 1 from Akametalu et al.,
# This version uses the k-nearest neighbors (KNN) method instead of interpolation.
# -----------------------------------------------------------------------------

# # --- Required Packages ---
# using Plots
# using Printf
# using LinearAlgebra
# using DelimitedFiles
# using NearestNeighbors # Add this for KNN functionality

# # ---------------------------
# # GLOBAL CONSTANTS
# # ---------------------------
# const X1_MIN, X1_MAX = -1.0, 5.0
# const X2_MIN, X2_MAX = -5.0, 5.0
# const NUM_POINTS_STATE_1 = 161
# const NUM_POINTS_STATE_2 = 161

# const U_MAX = 2.0
# const NUM_POINTS_ACTIONS = 41

# const K1_MIN, K1_MAX = 0.0, 4.0
# const K2_MIN, K2_MAX = -3.0, 3.0

# const LAMBDA = 0.2
# const TAU_BAR = 2.0
# const DT = 0.1
# const GAMMA = exp(-LAMBDA * DT)

# const TOLERANCE = 1e-6
# const MAX_ITER = 1000
# const L = sqrt((5.0 - 4.0)^2 + (5.0 - 3.0)^2)

# # ---------------------------
# # HELPER FUNCTIONS
# # ---------------------------
# function signed_distance_to_box(px, py)
#     dx = max(K1_MIN - px, 0.0, px - K1_MAX)
#     dy = max(K2_MIN - py, 0.0, py - K2_MAX)
#     dist_outside = sqrt(dx^2 + dy^2)
#     if dist_outside > 0.0
#         return dist_outside
#     else
#         dist_to_xmin = px - K1_MIN; dist_to_xmax = K1_MAX - px
#         dist_to_ymin = py - K2_MIN; dist_to_ymax = K2_MAX - py
#         return -min(dist_to_xmin, dist_to_xmax, dist_to_ymin, dist_to_ymax)
#     end
# end

# function dynamics(x1, x2, u)
#     x1_next = x1 + x2 * DT
#     x2_next = x2 + u * DT
#     return [x1_next, x2_next] # Return a vector for knn
# end

# # NOTE: The interpolate_bilinear function has been removed.

# # ---------------------------
# # MAIN SIMULATION FUNCTION
# # ---------------------------
# function solve_and_plot()
#     println("--- Starting MDR Value Iteration with KNN ---")

#     # --- Grid and Initialization ---
#     x1_grid = collect(LinRange(X1_MIN, X1_MAX, NUM_POINTS_STATE_1))
#     x2_grid = collect(LinRange(X2_MIN, X2_MAX, NUM_POINTS_STATE_2))
#     actions = collect(LinRange(-U_MAX, U_MAX, NUM_POINTS_ACTIONS))
    
#     # Create a list of all 2D state tuples for easy iteration
#     state_2d = [(x, v) for x in x1_grid for v in x2_grid]
#     nstates = length(state_2d)
#     println("Grid size: $NUM_POINTS_STATE_1 x $NUM_POINTS_STATE_2 = $nstates states.")
    
#     h_matrix = [-signed_distance_to_box(x, v) - L for x in x1_grid, v in x2_grid]
#     h_vec = vec(h_matrix)
#     U = copy(h_vec)
#     U_next = similar(U)

#     # --- NEW: Build KDTree for KNN lookups ---
#     # Convert the list of states into a matrix for the tree
#     states_matrix = hcat([collect(s) for s in state_2d]...)
#     tree = KDTree(states_matrix)
#     println("KDTree for KNN has been built.")

#     # --- Value Iteration Loop ---
#     println("\nStarting value iteration for λ = $LAMBDA...")
#     iteration = 0
#     diff = Inf
#     @time while diff > TOLERANCE && iteration < MAX_ITER
#         iteration += 1
#         max_diff = 0.0
        
#         # Note: We no longer need to reshape U to a matrix
        
#         for i in 1:nstates
#             x1, x2 = state_2d[i]
#             best_over_u = -Inf
            
#             for u in actions
#                 next_state = dynamics(x1, x2, u)
                
#                 # --- MODIFIED: Use KNN to find the value at the next state ---
#                 # Find the index of the single closest grid point (k=1)
#                 idxs, _ = knn(tree, next_state, 1)
#                 j = idxs[1]
                
#                 # Get the value of U at that closest grid point
#                 val_at_neighbor = GAMMA * U[j]
                
#                 best_over_u = max(best_over_u, val_at_neighbor)
#             end
            
#             U_next[i] = min(h_vec[i], best_over_u)
#             max_diff = max(max_diff, abs(U_next[i] - U[i]))
#         end
        
#         U .= U_next
#         diff = max_diff
        
#         if iteration % 10 == 0
#             @printf "Iteration %4d, max diff = %.8f\n" iteration diff
#         end
#     end

#     if diff <= TOLERANCE
#         println("Value iteration converged in $iteration iterations.")
#     else
#         @warn "Value iteration did not converge within $MAX_ITER iterations."
#     end
    
#     # --- Post-Processing and Saving ---
#     Z = reshape(U, (NUM_POINTS_STATE_1, NUM_POINTS_STATE_2)) .+ L
    
#     results_dir = joinpath(@__DIR__, "knn_results")
#     mkpath(results_dir)
#     output_path = joinpath(results_dir, "Value_function_lambda_0.2_knn.txt")
#     println("\nSaving Z matrix to: $output_path")
#     writedlm(output_path, Z, ',')
#     println("--- Data saved successfully. ---")
    
#     # --- Plotting ---
#     level_over = 0.0
#     level_under = L * (1 - exp(-LAMBDA * TAU_BAR))
#     @printf "\nPlotting Over-approx (Z(x)>%.2f) and Under-approx (Z(x)>%.4f)\n" level_over level_under

#     safe_box_shape = Shape([K1_MIN, K1_MAX, K1_MAX, K1_MIN], [K2_MIN, K2_MIN, K2_MAX, K2_MAX])
#     p = plot(safe_box_shape, opacity=0.3, color=:red, label="Safe Box K",
#         xlabel="State x₁", ylabel="State x₂",
#         title="MDR Safe Set Approximations (λ=$LAMBDA, KNN)",
#         aspect_ratio=:equal,
#         xlims=(X1_MIN, X1_MAX), ylims=(X2_MIN, X2_MAX),
#         legend=:bottomleft,
#         framestyle=:box
#     )
#     contour!(p, x1_grid, x2_grid, Z', levels=[level_over],
#         linewidth=3, color=:blue, linestyle=:solid, label="Z(x) Over-Approximation")
#     contour!(p, x1_grid, x2_grid, Z', levels=[level_under],
#         linewidth=3, color=:darkblue, linestyle=:dash, label="Z(x) Under-Approximation")

#     display(p)
#     println("--- Plot generated successfully. ---")
# end

# # ---------------------------
# # RUN THE SCRIPT
# # ---------------------------
# solve_and_plot()