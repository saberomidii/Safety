using Random, Distributions
using NearestNeighbors
using Plots
using DelimitedFiles # Added this as it was missing

# --- Setup: Automatically create a dummy policy file for this example ---
# This ensures the script is self-contained and runnable.
const NUM_STATES_TOTAL = 76 * 76 * 21
const NUM_ACTIONS = 21
# -----------------------------------------------------------------------

Random.seed!(10) # Setting the seed
cd(@__DIR__)
# Load the optimal policy.
policy_average = readdlm("optimal_policy_avr.csv", ';')
policy_MDR = readdlm("MDR_policy_map_lambda_0.0.csv", ';')

# --- Define Constants FIRST ---
# These must be defined before they are used in the functions below.
const a_outer = 0.4
const b_outer = 0.9
const a_inner = 0.2
const b_inner = 0.70
const x_c = 0.0
const v_c = 0.0 # This corresponds to the y-center for plotting

# --- Function Definitions ---
# Checks if the point (x, y) is within the safe ring.
function is_safe(x::Float64, y::Float64)
    inside_outer = ((x - x_c)^2) / a_outer^2 + ((y - v_c)^2) / b_outer^2 <= 1
    inside_inner = ((x - x_c)^2) / a_inner^2 + ((y - v_c)^2) / b_inner^2 <= 1
    return inside_outer && !inside_inner
end

# Continuous (noisy) dynamics for the Dubins car.
function di_dynamics(x1::Float64, x2::Float64, x3::Float64, u::Float64, d::Float64)
# if !is_safe(x1, x2)
#     return (x1, x2, x3) # No movement if outside the safe region
# else
    dt = 0.1
    V = 0.1 # Constant Speed
    x1_next = x1 + V * cos(x3) * dt
    # CORRECTED: The y-component (x2) must use sin(theta).
    x2_next = x2 + V * sin(x3) * dt
    x3_next = x3 + (u + d) * dt
    return (x1_next, x2_next, x3_next)
end

# --- State Space Discretization ---
const num_points_state_1 = 76
const num_points_state_2 = 76
const num_points_state_3 = 21
const x1_range = collect(LinRange(-0.5, 0.5, num_points_state_1))
const x2_range = collect(LinRange(-0.5, 0.5, num_points_state_2))
const x3_range = collect(LinRange(0, 2π, num_points_state_3))

const states_3d = [(x, v, theta) for x in x1_range for v in x2_range for theta in x3_range]
const nstates = length(states_3d)
const states_matrix = hcat([collect(s) for s in states_3d]...)
const tree = KDTree(states_matrix)

# --- Simulation ---
const sim_time = 200 # Increased for a more visible trajectory
const μ = 0.0
const σ = 1.0
const low = -1.0
const up = 1.0
const u = collect(LinRange(-1.0, 1.0, NUM_ACTIONS))

# Initial state
xn = 0.25
yn = 0.25 # Renamed from 'vn' for clarity in X-Y plot
thetan = 0.5

# Store history as vectors
X = zeros(sim_time)
Y = zeros(sim_time)
Theta = zeros(sim_time)

println("Running simulation...")
for time_index = 1:sim_time
    # Find the nearest discrete state
    idxs, _ = knn(tree, [xn, yn, thetan], 1)
    
    # Get the policy for that state
    policy_index = policy_average[idxs[1]]
    policy_index = policy_MDR[idxs[1]]

    policy_index = round(Int64, policy_index)
    
    # Generate a random disturbance
    disturbance = rand(Truncated(Normal(μ, σ), low, up))
    
    # FIXED: Use `global` to modify variables outside the loop's scope
    global (xn, yn, thetan) = di_dynamics(xn, yn, thetan, u[policy_index], disturbance)
    
    # Store the results
    X[time_index] = xn
    Y[time_index] = yn
    Theta[time_index] = thetan
end
println("Simulation finished.")

# --- Plotting ---
# Helper function to generate points for an ellipse
function generate_ellipse_points(a, b, xc=0.0, yc=0.0; n_points=100)
    t = range(0, 2π, length=n_points)
    x_points = xc .+ a .* cos.(t)
    y_points = yc .+ b .* sin.(t)
    return x_points, y_points
end

# Generate points for both ellipses
outer_ellipse_x, outer_ellipse_y = generate_ellipse_points(a_outer, b_outer, x_c, v_c)
inner_ellipse_x, inner_ellipse_y = generate_ellipse_points(a_inner, b_inner, x_c, v_c)

# Create the plot
plot(X, Y,
    label="Dubins Car Trajectory",
    xlabel="x position",
    ylabel="y position",
    title="Dubins Car Simulation with Safe Set",
    legend=:topright,
    linewidth=2,
    aspect_ratio=:equal, # Ensures ellipses are not distorted
    framestyle=:box
)

# Add the outer ellipse and shade the safe area
plot!(outer_ellipse_x, outer_ellipse_y,
    label="Safe Set Boundary",
    seriestype=:shape, # Use shape to fill
    fillcolor=:green,
    fillalpha=0.2, # Transparent fill
    linecolor=:gray,
    linestyle=:dash
)

# Add the inner ellipse and "punch a hole" in the fill
plot!(inner_ellipse_x, inner_ellipse_y,
    label="Unsafe Inner Region",
    seriestype=:shape,
    fillcolor=:white, # Fill with background color
    fillalpha=1.0,
    linecolor=:red
)

# Mark the start point
scatter!([X[1]], [Y[1]], label="Start Point", color=:blue, markersize=5)

println("Displaying plot...")
display(current()) # Show the plot
