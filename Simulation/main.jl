# This script runs a robustness test for the double integrator system.
# It loads a pre-computed binary policy from the "Average Reward" results.

using Plots
using Distributions
using Random
using NearestNeighbors
using DelimitedFiles
using JLD2 # Required for reading the binary .jld2 file

# ---------------------------
# SIMULATION PARAMETERS
# ---------------------------
const NUM_SIMULATIONS = 50       # N: The number of runs to attempt
const DT = 0.1                   # Time step in seconds
const SIMULATION_STEPS = 150     # Max steps per run
const INITIAL_STATE = [3.0, 0.0] # Start at position 0, velocity 0

# --- Constraint Box K ---
const K1_MIN, K1_MAX = 0.0, 4.0
const K2_MIN, K2_MAX = -3.0, 3.0

# ---------------------------
# POLICY GRID PARAMETERS (must match the file)
# ---------------------------
const X1_MIN, X1_MAX = -1.0, 5.0
const X2_MIN, X2_MAX = -5.0, 5.0
const NUM_POINTS_STATE_1 = 161
const NUM_POINTS_STATE_2 = 161 # This must match the dimensions of the policy file

# --- ACTION SPACE (must match the file generation script) ---
const U_MIN, U_MAX = -2.0, 2.0
const NUM_POINTS_ACTIONS = 41

# ---------------------------
# DISTURBANCE PARAMETERS
# ---------------------------
const SIGMA = 1.0
const MEAN = 0.0
const D_MIN = -1.0
const D_MAX = 1.0

# ---------------------------
# HELPER FUNCTIONS
# ---------------------------
function is_safe(state)
    x1, x2 = state
    return (K1_MIN <= x1 <= K1_MAX) && (K2_MIN <= x2 <= K2_MAX)
end

function dynamics(current_state, u, d, dt)
    x1, x2 = current_state
    x1_next = x1 + dt * x2
    x2_next = x2 + dt * (u + d)
    return [x1_next, x2_next]
end

function optimal_policy(state_index, policy_vector, actions)
    action_index = policy_vector[state_index]
    return actions[action_index]
end

# ---------------------------
# MAIN SIMULATION FUNCTION
# ---------------------------
function run_robustness_test()
    println("--- Starting robustness test with loaded Average Reward policy ---")

    # --- ACTION SPACE AND POLICY LOADING ---
    actions = LinRange(U_MIN, U_MAX, NUM_POINTS_ACTIONS)
    
    # !!! IMPORTANT: UPDATE THIS PATH to your .jld2 file's location !!!
    policy_path = joinpath(@__DIR__, "results_primal", string(SIGMA), "optimal_policy.jld2")
    
    println("Loading binary policy from: $policy_path")
    
    # --- MODIFIED: Load the binary .jld2 file ---
    policy_data = jldload(policy_path)
    policy_map = policy_data["policy_map"] # Extracts the variable named "policy_map"
    policy_vector = vec(policy_map')       # Convert to vector, transposing to match Julia's order
    
    if length(policy_vector) != (NUM_POINTS_STATE_1 * NUM_POINTS_STATE_2)
        @warn "Dimension mismatch! Policy file size does not match grid size."
        return
    end

    # --- Create the policy grid and KDTree ---
    println("Building policy grid for KNN...")
    x1_grid = LinRange(X1_MIN, X1_MAX, NUM_POINTS_STATE_1)
    x2_grid = LinRange(X2_MIN, X2_MAX, NUM_POINTS_STATE_2)
    state_grid = [[x, v] for x in x1_grid for v in x2_grid]
    states_matrix = hcat(state_grid...)
    tree = KDTree(states_matrix)
    println("KDTree built successfully.")
    
    all_trajectories = []
    simulation_outcomes = []
    first_failure_index = 0

    # --- Main loop over N simulations ---
    for i in 1:NUM_SIMULATIONS
        current_trajectory = [INITIAL_STATE]
        current_state = INITIAL_STATE
        run_failed = false

        for k in 1:SIMULATION_STEPS
            idxs, _ = knn(tree, current_state, 1)
            nearest_index = idxs[1]
            u = optimal_policy(nearest_index, policy_vector, actions)
            d = clamp(rand(Normal(MEAN, SIGMA)), D_MIN, D_MAX)
            next_state = dynamics(current_state, u, d, DT)
            push!(current_trajectory, next_state)
            current_state = next_state
            if !is_safe(current_state)
                run_failed = true
                break
            end
        end
        
        push!(all_trajectories, current_trajectory)
        if run_failed
            push!(simulation_outcomes, :failed)
            first_failure_index = i
            println("\n!!! Simulation FAILED on run #$i. Stopping test. !!!")
            break
        else
            push!(simulation_outcomes, :safe)
            println("Run #$i completed safely.")
        end
    end

    println("\n--- Test finished. Plotting results... ---")
    
    # --- PLOTTING LOGIC ---
    box_shape = Shape([K1_MIN, K1_MAX, K1_MAX, K1_MIN], [K2_MIN, K2_MIN, K2_MAX, K2_MAX])
    
    if first_failure_index > 0
        p = plot(box_shape, opacity=0.3, color=:gray, label="Constraint Box",
            title="First Failure on Run #$first_failure_index / $NUM_SIMULATIONS",
            xlabel="State x₁ (Position)", ylabel="State x₂ (Velocity)")
        for (idx, trajectory) in enumerate(all_trajectories)
            x1_traj = [s[1] for s in trajectory]; x2_traj = [s[2] for s in trajectory]
            plot!(p, x1_traj, x2_traj, color=(simulation_outcomes[idx] == :safe ? :green : :red), label="")
        end
        plot!(p, [], [], color=:green, label="Safe Runs"); plot!(p, [], [], color=:red, label="Failed Run")
    else
        p = plot(box_shape, opacity=0.3, color=:gray, label="Constraint Box",
            title="$NUM_SIMULATIONS Successful Runs - Showing Last Trajectory",
            xlabel="State x₁ (Position)", ylabel="State x₂ (Velocity)")
        last_trajectory = all_trajectories[end]
        x1_traj = [s[1] for s in last_trajectory]; x2_traj = [s[2] for s in last_trajectory]
        plot!(p, x1_traj, x2_traj, color=:green, label="Successful Run", marker=:circle)
    end
    
    display(p)
end

# --- RUN THE SCRIPT ---
run_robustness_test()