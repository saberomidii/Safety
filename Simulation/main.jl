using Random
using Distributions
using NearestNeighbors
using DelimitedFiles

# --- Setup: Ensure a dummy policy file exists ---
function create_dummy_policy_file(filename, num_states, num_actions)
    if !isfile(filename)
        println("Creating dummy policy file: $filename")
        dummy_policy = [round(Int, num_actions / 2) for _ in 1:num_states]
        writedlm(filename, dummy_policy, ';')
    end
end

const NUM_STATES_TOTAL = 61 * 61 * 31
const NUM_ACTIONS = 41
create_dummy_policy_file("optimal_policy_avr.csv", NUM_STATES_TOTAL, NUM_ACTIONS)
# -----------------------------------------------------------------------

Random.seed!(10) # For reproducible results
cd(@__DIR__) # Set working directory to the script's location

println("Loading policy and setting up state space...")
const policy_average = readdlm("optimal_policy_avr.csv", ';')

# --- ENVIRONMENT AND DYNAMICS DEFINITIONS ---
const L_straight = 2.0
const R_outer = 1.75
const R_inner = 1.0
const xc_left = -L_straight / 2.0
const xc_right = L_straight / 2.0
const y_c = 0.0

# --- HELPER FUNCTIONS ---
function is_point_safe(x::Float64, y::Float64)
    is_in_straights = (x >= xc_left && x <= xc_right) && (abs(y - y_c) <= R_outer)
    is_in_left_curve = x < xc_left && ((x - xc_left)^2 + (y - y_c)^2 <= R_outer^2)
    is_in_right_curve = x > xc_right && ((x - xc_right)^2 + (y - y_c)^2 <= R_outer^2)
    inside_outer = is_in_straights || is_in_left_curve || is_in_right_curve
    
    is_in_inner_straights = (x >= xc_left && x <= xc_right) && (abs(y - y_c) < R_inner)
    is_in_inner_left_curve = x < xc_left && ((x - xc_left)^2 + (y - y_c)^2 < R_inner^2)
    is_in_inner_right_curve = x > xc_right && ((x - xc_right)^2 + (y - y_c)^2 < R_inner^2)
    inside_inner = is_in_inner_straights || is_in_inner_left_curve || is_in_inner_right_curve

    return inside_outer && !inside_inner
end

function di_dynamics(x1, x2, x3, u, d)
    dt = 0.1
    V = .2
    x1_next = x1 + V * cos(x3) * dt
    x2_next = x2 + V * sin(x3) * dt
    x3_next = x3 + (u + d) * dt
    return (x1_next, x2_next, x3_next)
end

function generate_racetrack_points(R, L; y_c=0.0, n_points_curve=100)
    if R <= 0 return Float64[], Float64[] end
    xc_l, xc_r = -L / 2.0, L / 2.0
    t_right = LinRange(π/2, -π/2, n_points_curve)
    x_r_curve = xc_r .+ R .* cos.(t_right)
    y_r_curve = y_c .+ R .* sin.(t_right)
    t_left = LinRange(3π/2, π/2, n_points_curve)
    x_l_curve = xc_l .+ R .* cos.(t_left)
    y_l_curve = y_c .+ R .* sin.(t_left)
    # Combine to make a closed shape for filling
    x_points = vcat(x_r_curve, xc_l, x_l_curve, xc_r)
    y_points = vcat(y_r_curve, -R, y_l_curve, R)
    return x_points, y_points
end

# --- STATE SPACE DISCRETIZATION ---
println("Discretizing state space and building KD-Tree...")
const x1_range = LinRange(-3.0, 3.0, 61)
const x2_range = LinRange(-3.0, 3.0, 61)
const x3_range = LinRange(0, 2π, 31)
const states_3d = [(x, y, θ) for x in x1_range for y in x2_range for θ in x3_range]
const states_matrix = hcat(collect.(states_3d)...)
const tree = KDTree(states_matrix)

# --- SIMULATION PARAMETERS ---
const u_actions = LinRange(-1.75, 1.75, NUM_ACTIONS)
const max_sim_time = 1000
const disturbance_dist = Truncated(Normal(0.0, 1.0), -1., 1.) # Define distribution once


# --- SIMULATION FUNCTION ---
function run_simulation(x_start, y_start, theta_start)
    X_traj, Y_traj, U_traj = Float64[], Float64[], Float64[]
    xn, yn, thetan = x_start, y_start, theta_start
    
    for _ in 1:max_sim_time
        is_point_safe(xn, yn) || break
        push!(X_traj, xn); push!(Y_traj, yn)
        
        idxs, _ = knn(tree, [xn, yn, thetan], 1)
        policy_index = round(Int, policy_average[idxs[1]])
        input_u = u_actions[policy_index]
        push!(U_traj, input_u)

        disturbance = rand(disturbance_dist)
        (xn, yn, thetan) = di_dynamics(xn, yn, thetan, input_u, disturbance)
    end
    
    return X_traj, Y_traj, U_traj
end

# # --- MAIN EXECUTION ---
# longest_time = 0
# best_starting_point = (0.0, 0.0, 0.0)
# total_states_to_check = length(states_3d)

# println("\n--- Searching for the best starting point out of $total_states_to_check candidates ---")

# for (i, start_state) in enumerate(states_3d)
#     is_point_safe(start_state[1], start_state[2]) || continue
#     if i % 20000 == 0
#         println("... Progress: $(round(100 * i / total_states_to_check, digits=1))%")
#     end

#     X_traj, _, _ = run_simulation(start_state...)
#     current_time_safe = length(X_traj)

#     if current_time_safe > longest_time
#         longest_time
#         global longest_time = current_time_safe
#         global best_starting_point = start_state
#         println("✅ New best found! Time: $longest_time steps, Start: (x=$(round(start_state[1], digits=2)), y=$(round(start_state[2], digits=2)), θ=$(round(start_state[3], digits=2)))")
#     end
# end

# println("\n--- 🏆 Search Complete! ---")
# println("Best Starting Point: $best_starting_point")
# println("Longest Time Safe: $longest_time steps")

# --- SAVE DATA FOR PLOTTING ---
println("\n--- Saving data for MATLAB plotting ---")
data_dir = "figure_data"
mkpath(data_dir) # Create directory if it doesn't exist

# Run final simulation to get the best trajectory
(best_x, best_y, best_theta) = (0.0, -1.5, pi)

# (best_x, best_y, best_theta) = best_starting_point

best_X_traj, best_Y_traj, best_U_traj = run_simulation(best_x, best_y, best_theta) 

# Generate racetrack points
outer_track_x, outer_track_y = generate_racetrack_points(R_outer, L_straight)
inner_track_x, inner_track_y = generate_racetrack_points(R_inner, L_straight)

# Save data files
writedlm(joinpath(data_dir, "trajectory_3.csv"), hcat(best_X_traj, best_Y_traj), ',')
writedlm(joinpath(data_dir, "control_input_3.csv"), best_U_traj, ',')
writedlm(joinpath(data_dir, "outer_track.csv"), hcat(outer_track_x, outer_track_y), ',')
writedlm(joinpath(data_dir, "inner_track.csv"), hcat(inner_track_x, inner_track_y), ',')
writedlm(joinpath(data_dir, "start_point_3.csv"), [best_x, best_y], ',')
# writedlm(joinpath(data_dir, "longest_time.txt"), longest_time)

println("Data saved to '$data_dir' directory.")