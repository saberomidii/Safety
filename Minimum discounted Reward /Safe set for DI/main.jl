using Plots
using Printf
using LinearAlgebra
using DelimitedFiles
using NearestNeighbors
using Distributions
using SparseArrays
using Random
using DataFrames

Random.seed!(10)



# Above Left Corner. 1
# const x_1_min, x_1_max = -0.5, 1.5
# const x_2_min, x_2_max = 1.5, 3.5

## Bottom Left Corner
# const x_1_min, x_1_max = -0.5, 1.5
# const x_2_min, x_2_max = -3.5, -1.5


# # Bottom Right Corner. 2
# const x_1_min, x_1_max = 2.5, 4.5
# const x_2_min, x_2_max = -3.5, -1.5

## Above Right Corner
# const x_1_min, x_1_max = 2.5, 4.5
# const x_2_min, x_2_max = 1.5, 3.5


const x_1_min, x_1_max = -1.0, 5.0
const x_2_min, x_2_max = -5.0, 5.0

const u_min, u_max = -2.0, 2.0


const num_points_state_1 = 2*160+1
const num_points_state_2 = 2*160+1

const num_points_action = 81


const K1_MIN, K1_MAX = 0.0, 4.0
const K2_MIN, K2_MAX = -3.0, 3.0

const nsamples = 100
const sigma = 1.0
const threshold_for_transit = 0.000

### Read disturbance
script_dir = @__DIR__
safety_folder_path = joinpath(script_dir, "..", "..")
csv_filepath = joinpath(safety_folder_path, "Disturbance.csv")
disturbance_list = readdlm(csv_filepath, ',', Float64)

# Number of states in each dimension
x1 = collect(LinRange(x_1_min, x_1_max, num_points_state_1))  # Possible x values
x2 = collect(LinRange(x_2_min, x_2_max, num_points_state_2))  # Possible v values

# Cartesian product of (x, v) forms the entire state space
const states_2d = [(x, v) for x in x1 for v in x2]
const nstates = length(states_2d)

# Action discretization
const actions = collect(LinRange(u_min, u_max, num_points_action))
const nactions = length(actions)

function is_safe(x::Float64, v::Float64)
    # "Safe" region: x in [0,4], v in [-3,3]
    return ( K1_MIN<= x <= K1_MAX) && (K2_MIN <= v <= K2_MAX)
end

function signed_distance_to_box(px, py)
    dx = max(K1_MIN - px, 0.0, px - K1_MAX)
    dy = max(K2_MIN - py, 0.0, py - K2_MAX)
    dist_outside = sqrt(dx^2 + dy^2) 
    return dist_outside > 0 ? dist_outside : -min(px - K1_MIN, K1_MAX - px, py - K2_MIN, K2_MAX - py)
end

# Continuous (noisy) dynamics for double-integrator
function dynamics_rand(x::Float64, v::Float64, u::Float64, d::Float64)    
    if !is_safe(x, v)
        return (x, v)   # no movement if outside the safe region
    else
        dt=0.1
        x1_next = x + v*dt
	    x2_next = v + (u + d)*dt
        return (x1_next, x2_next)
    end
end


# Build a KD-tree for snapping continuous next-states to the nearest discrete state
states_matrix = hcat([collect(s) for s in states_2d]...)  # shape: 2 x nstates
tree = KDTree(states_matrix)

println("\nBuilding transition probabilities T[s][a, s_next] ...")

# Initialize T array: T[s][a, s_next]
T= Vector{SparseMatrixCSC{Float64,Int64}}(undef,nstates)


# Initialize T array: T[s][a, s_next]

for index_state in 1:nstates
    T[index_state] = spzeros(nactions,nstates)
end

# (Compute Transition Matrix)
for is in 1:nstates
    s = states_2d[is]   # e.g., s = (x, v)
    for a in 1:nactions
        for i in 1:nsamples
            xn = s[1]
            vn = s[2]
            j = 1
            while true	
                random_index=rand(1:nsamples)
                (xn, vn) = dynamics_rand(s[1], s[2], actions[a],disturbance_list[random_index])
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

println("Negative values checking")
@assert all(t->all(t.>=0.0),T) "Negative Value!"

println("Rows sums are approximately 1!")
for action in 1:num_points_action
    @assert all(i -> sum(T[i][action, :]) ≈ 1.0, axes(T, 1)) "Not All row sums are approximately1!"
end

println("Removing small values")
for s in 1:nstates
    T[s]= dropzeros!(map(x->abs(x)<threshold_for_transit ? 0.0 : x, T[s]))
end

println("No Error")

#_____ MDR section ________
const LAMBDA = 0.0
const TAU_BAR = 2.0
const DT = 0.1
const time_step = 1000
const TOLERANCE = 1e-9 # A small number to check against
converged_t = 1 # Default value if the loop completes



const GAMMA = exp(-LAMBDA * DT)
const TOLERANCE = 1e-16
const MAX_ITER = 10000
const L = sqrt((5.0 - 4.0)^2 + (5.0 - 3.0)^2)

function signed_distance_to_box(px, py)
    dx = max(K1_MIN - px, 0.0, px - K1_MAX)
    dy = max(K2_MIN - py, 0.0, py - K2_MAX)
    dist_outside = sqrt(dx^2 + dy^2) 
    return dist_outside > 0 ? dist_outside : -min(px - K1_MIN, K1_MAX - px, py - K2_MIN, K2_MAX - py)
end



x1_grid = LinRange(x_1_min, x_1_max, num_points_state_1)
x2_grid = LinRange(x_2_min, x_2_max, num_points_state_2)
h_matrix = [-signed_distance_to_box(x, v) - L for x in x1_grid, v in x2_grid]


println("--- Step 2: Solving MDR using 2D matrices for states and values  for simple case---")
# Value function is kept as a 2D matrix

U = zeros(Float64, time_step, num_points_state_1, num_points_state_2)

U[end,:,:] = h_matrix

for t in time_step-1:-1:1    
    # Loop explicitly over the 2D grid indices
    for iy in 1:num_points_state_2
        for ix in 1:num_points_state_1
            best_over_u = -Inf

            # Convert 2D index (ix, iy) to 1D linear index to access T
            # linear_idx_current = (iy - 1) * num_points_state_1 + ix
            linear_idx_current = (ix - 1) * num_points_state_2 + iy



            for a_idx in 1:num_points_action
                # Find all possible next states (as 1D indices) from T
                s_primes_indices_1D = findnz(T[linear_idx_current][a_idx, :])[1]
                
                if isempty(s_primes_indices_1D); continue; end
                
                min_val_of_next_state = Inf
                # For each possible next state, find its value in the 2D U_matrix
                for next_state_1D_idx in s_primes_indices_1D
                    # Convert the 1D next_state index back to 2D indices
                    # ix_prime = (next_state_1D_idx - 1) % num_points_state_1 + 1
                    ix_prime = fld(next_state_1D_idx - 1, num_points_state_2) + 1

                    # iy_prime = fld(next_state_1D_idx - 1, num_points_state_1) + 1
                    iy_prime = (next_state_1D_idx - 1) % num_points_state_2 + 1

                    val = U[t+1,ix_prime, iy_prime]
                    min_val_of_next_state = min(min_val_of_next_state, val)
                end
                
                worst_case_value = min_val_of_next_state
                best_over_u = max(best_over_u, GAMMA * worst_case_value)
            end
            
            U[t, ix, iy] = min(h_matrix[ix, iy], best_over_u)
        end
    end
    # --- ADD THE CONVERGENCE CHECK CODE HERE ---
        diff = maximum(abs.(U[t,:,:] .- U[t+1,:,:]))
        println("Time step: ", t, ", Max change: ", diff)

        if diff < TOLERANCE
            println("\nValue function converged at time step t=", t)
            converged_t = t # Save the current time step

            for i in 1:t-1
                U[i,:,:] = U[t,:,:]
            end
            break
        end
end

# Calculate the final Z value from the converged U matrix
Z = (U .+ L)

cd(@__DIR__)
writedlm("Z_function.csv", Z[converged_t,:,:]', ',')

# println("--- Generating Plot ---")

# # --- Dynamic Plotting Logic for Any View ---

# # Create the base contour plot for the zero level-set.
# # I've enabled the legend to make sure labels appear.
# p = contour(x1_grid, x2_grid, Z[converged_t,:,:]',
#             levels=[0],
#             color=:blue,
#             lw=2.5,
#             label="Zero Level-Set",
#             legend=:topright)

# # Initialize variables for the plot title and safe set coordinates
# plot_title = "State Space View"
# safe_x_coords = []
# safe_y_coords = []

# # Define the absolute boundaries of the safe set
# const K1_MIN, K1_MAX = 0.0, 4.0
# const K2_MIN, K2_MAX = -3.0, 3.0

# # Logic to handle different views
# if x_1_min == -0.5 && x_2_min == 1.5
#     # Case: Top Left Corner
#     plot_title = "Top Left Corner"
#     safe_x_coords = [0.0, 1.5, 1.5, 0.0]
#     safe_y_coords = [1.5, 1.5, 3.0, 3.0]

# elseif x_1_min == 2.5 && x_2_min == -3.5
#     # Case: Bottom Right Corner
#     plot_title = "Bottom Right Corner"
#     safe_x_coords = [2.5, 4.0, 4.0, 2.5]
#     safe_y_coords = [-3.0, -3.0, -1.5, -1.5]

# # Add other specific corner cases here if you need them...

# else
#     # General Case: This will handle your new, larger view [-1,5]x[-5,5]
#     plot_title = "Full State Space View"
#     # For a general view, we draw the entire safe set rectangle
#     safe_x_coords = [K1_MIN, K1_MAX, K1_MAX, K1_MIN] # [0, 4, 4, 0]
#     safe_y_coords = [K2_MIN, K2_MIN, K2_MAX, K2_MAX] # [-3, -3, 3, 3]
# end

# # Add the safe set shape to the plot.
# # The light green shaded area indicates the safe set itself.
# if !isempty(safe_x_coords)
#     plot!(p, safe_x_coords, safe_y_coords,
#           seriestype=:shape,
#           c=:lightgreen,
#           alpha=0.4,
#           linecolor=nothing,
#           label="Safe Set")
# end

# # Add the final title and axis labels
# title!(p, plot_title)
# xlabel!(p, "x₁")
# ylabel!(p, "x₂")

# # Display the final plot
# display(p)
