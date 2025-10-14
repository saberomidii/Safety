using Random
using Distributions
using NearestNeighbors
using Plots
using DelimitedFiles

# --- Setup: Ensure a dummy policy file exists ---
function create_dummy_policy_file(filename, num_states, num_actions)
    if !isfile(filename)
        println("Creating dummy policy file: $filename")
        dummy_policy = [round(Int, num_actions / 2) for _ in 1:num_states]
        writedlm(filename, dummy_policy, ';')
    end
end

const NUM_STATES_TOTAL = 81 * 81 * 21
const NUM_ACTIONS = 81
create_dummy_policy_file("optimal_policy_avr.csv", NUM_STATES_TOTAL, NUM_ACTIONS)
# -----------------------------------------------------------------------

Random.seed!(10) # For reproducible resuls
cd(@__DIR__) # Set working directory to the script's location

println("Loading policy and setting up state space...")
const policy_average = readdlm("optimal_policy_avr.csv", ';')

# --- Function Definitions ---
# --- Define Racetrack Constants ---
const L_straight = 0.9 # Length of the straight sections
const R_outer = 0.4    # Radius of the outer semicircles
const R_inner = 0.2    # Radius of the inner semicircles
const xc_left = -L_straight / 2.0  # x-center of the left semicircle
const xc_right = L_straight / 2.0 # x-center of the right semicircle
const y_c = 0.0 # y-center is at 0

# --- Function Definitions ---

# Helper function to check if a point is inside a single racetrack shape
function is_inside_racetrack(x::Float64, y::Float64, R::Float64)
    # Check straight sections
    if x >= xc_left && x <= xc_right
        return abs(y - y_c) <= R
    # Check left semicircle
    elseif x < xc_left
        return (x - xc_left)^2 + (y - y_c)^2 <= R^2
    # Check right semicircle
    else # x > xc_right
        return (x - xc_right)^2 + (y - y_c)^2 <= R^2
    end
end

# Checks if a point (x, y) is in the safe region between two racetrack shapes
function is_point_safe(x::Float64, y::Float64)
    inside_outer = is_inside_racetrack(x, y, R_outer)
    inside_inner = is_inside_racetrack(x, y, R_inner)
    return inside_outer && !inside_inner
end

# Calculates the next state based on dynamics
function di_dynamics(x1::Float64, x2::Float64, x3::Float64, u::Float64, d::Float64)
    dt = 0.1
    V =1.0 # Constant Speed
    x1_next = x1 + V * cos(x3) * dt
    x2_next = x2 + V * sin(x3) * dt
    x3_next = x3 + (u + d) * dt
    return (x1_next, x2_next, x3_next)
end

# --- State Space Discretization ---
println("Discretizing state space and building KD-Tree...")
const x1_range = LinRange(-1.0, 1.0, 81)
const x2_range = LinRange(-0.5, 0.5, 81)
const x3_range = LinRange(0, 2π, 21)

const states_3d = [(x, y, θ) for x in x1_range for y in x2_range for θ in x3_range]
const states_matrix = hcat(collect.(states_3d)...)
const tree = KDTree(states_matrix)

# --- Simulation Parameters ---
const μ = 0.0
const σ = 1.0
const disturbance_distribution = Truncated(Normal(μ, σ), -1.0, 1.0)
const u_actions = LinRange(-2.0, 2.0, NUM_ACTIONS)
const max_sim_time = 40# Set a reasonable maximum to prevent infinite loops

# --- Simplified Simulation Function ---
function run_simulation(x_start::Float64, y_start::Float64, theta_start::Float64)
    # Use dynamically sized vectors for the trajectory
    X_traj, Y_traj = Float64[], Float64[]
    xn, yn, thetan = x_start, y_start, theta_start
    
    # Loop until the agent is unsafe or max time is reached
    for _ in 1:max_sim_time
        # !is_point_safe(xn, yn) && break # Exit if unsafe

        push!(X_traj, xn)
        push!(Y_traj, yn)
        
        # Find the closest state in our grid to get the policy action
        idxs, _ = knn(tree, [xn, yn, thetan], 1)
        policy_index = round(Int, policy_average[idxs[1]])
        
        # Apply action with random disturbance
        disturbance = rand(disturbance_distribution)
        (xn, yn, thetan) = di_dynamics(xn, yn, thetan, u_actions[policy_index], disturbance)
    end
    
    return X_traj, Y_traj
end

# --- Main Execution ---
println("\n--- Running a single simulation from a random initial state ---")


p = plot(
    title="Agent Trajectory",
    xlabel="x position",
    ylabel="y position",
    legend=:outertopright,
    aspect_ratio=:equal,
    framestyle=:box
)

# --- Plotting ---
function generate_racetrack_points(R, L; y_c=0.0, n_points_curve=50)
    xc_l, xc_r = -L / 2.0, L / 2.0
    
    # Points for right semicircle
    t_right = LinRange(π/2, -π/2, n_points_curve)
    x_r_curve = xc_r .+ R .* cos.(t_right)
    y_r_curve = y_c .+ R .* sin.(t_right)
    
    # Points for left semicircle
    t_left = LinRange(3π/2, π/2, n_points_curve)
    x_l_curve = xc_l .+ R .* cos.(t_left)
    y_l_curve = y_c .+ R .* sin.(t_left)
    
    # Combine points to trace the full perimeter
    return vcat(x_r_curve, x_l_curve), vcat(y_r_curve, y_l_curve)
end

outer_track_x, outer_track_y = generate_racetrack_points(R_outer, L_straight)
inner_track_x, inner_track_y = generate_racetrack_points(R_inner, L_straight)

# Plot safe and unsafe regions
plot!(p, outer_track_x, outer_track_y, seriestype=:shape, c=:gray, alpha=0.2, label="Safe Region")
plot!(p, inner_track_x, inner_track_y, seriestype=:shape, c=:white, alpha=1.0, lc=:red, label="Unsafe Region")



# Run simulations from multiple safe starting points
for x_index in x1_range
    for y_index in x2_range
        for theta_index in x3_range
            if is_point_safe(x_index, y_index)
                X_traj, Y_traj = run_simulation(x_index, y_index, theta_index)

                # Check if the entire trajectory remained within the safe zone
                if all(is_point_safe.(X_traj, Y_traj))
                    println("Found a safe trajectory starting at: ", [x_index, y_index, theta_index])
                    plot!(p, X_traj, Y_traj, lw=2, color=:blue, legend=false)
                    # We can break after finding one to avoid a cluttered plot, or let it run
                end
            end
        end
    end
end

display(p)


