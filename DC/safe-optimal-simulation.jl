# Safe Optimal simulation for dubins path. 
# Written by Saber Omidi 

# Import required libraries
using Random
using Distributions
using NearestNeighbors
using SparseArrays
using LinearAlgebra
using JuMP
using Mosek
using MosekTools
using Plots
using Dates
using DelimitedFiles

# Constant Seed for Randomness 
Random.seed!(42)

# Grid Size  a rectangular x,y and theta for orientation
num_points_x = 55
num_points_y = 55
num_points_th = 32

## x y range    
# Rectangular area: 6m long (x) by 4m wide (y)
x_grid = collect(LinRange(-3.0, 3.0, num_points_x))
y_grid = collect(LinRange(-3.0, 3.0, num_points_y))
th_grid = collect(LinRange(0, 2π, num_points_th))

# Race track as a constrained set.
function constraint_set(x_grid, y_grid)
    cx = Float64[]
    cy = Float64[]
    
    for x in x_grid
        for y in y_grid
            # Geometry Logic: Distance from the central "skeleton" line
            dx_rect = max(abs(x) - L_straight/2.0, 0.0)
            dy_rect = abs(y - y_c)
            d_center = sqrt(dx_rect^2 + dy_rect^2)

            # If inside the track boundaries, add to the set
            if (d_center <= R_outer) && (d_center >= R_inner)
                push!(cx, x)
                push!(cy, y)
            end
        end
    end
    
    return cx, cy
end

L_straight = 2.0 
R_outer    = 1.80    
R_inner    = 0.90    
xc_left    = -L_straight / 2.0
xc_right   = L_straight / 2.0
y_c        = 0.0

## plot to test the constraint set 
# 1. Get the set of valid points
cx, cy = constraint_set(x_grid, y_grid)

# # 2. Plot
# p=scatter(cx, cy, 
#     aspect_ratio=:equal, 
#     markersize=1.5, 
#     color=:black,
#     legend=false, 
#     title="Constraint Set: Race Track"
# )

# Adding obstacles in the race track to make risky descions 
# === OBSTACLE ===
Obs_x = 0.0
Obs_y = 1.1 # Fixed location
Obs_r = 0.4

# === OBSTACLE FUNCTION ===
function obstacle_set(x_grid, y_grid)
    ox, oy = Float64[], Float64[]
    for x in x_grid, y in y_grid
        if (x - Obs_x)^2 + (y - Obs_y)^2 <= Obs_r^2
            push!(ox, x); push!(oy, y)
        end
    end
    return ox, oy
end

ox, oy = obstacle_set(x_grid, y_grid)

# Ploting to check to Obstacle location 
#p = scatter(cx, cy, c=:black, ms=1.5, aspect_ratio=:equal, label="Track", title="Track & Obstacle")

# # 3. Scatter Plot Obstacle (Red)
# scatter!(p, ox, oy, c=:red, ms=1.5, label="Obstacle")
# display(p)


# Is safe function returns bolean one or zero for those are inside the constraint set and they are not inside the obstable. 
function is_safe(x, y)
    # 1. Track Boundary Check
    dx = max(abs(x) - 1.0, 0.0); dy = abs(y)
    d = sqrt(dx^2 + dy^2)
    in_track = (d <= R_outer) && (d >= R_inner)
    
    # 2. Obstacle Check
    in_obs = (x - Obs_x)^2 + (y - Obs_y)^2 <= Obs_r^2
    
    return in_track && !in_obs
end 

# ploting to check zero and one  bolean 
heatmap_data = [is_safe(x, y) ? 1.0 : 0.0 for y in y_grid, x in x_grid]

p = heatmap(x_grid, y_grid, heatmap_data, 
    c=:grays, 
    aspect_ratio=:equal, 
    title="Safety Map: 1=Safe (White), 0=Unsafe (Black)",
    xlabel="X (m)", ylabel="Y (m)"
)
START_STATE  = [-2.5, -0.15, pi/16] 
TARGET_STATE = [1.5, 1.35, 0] 

# --- 4.2 Adding Start and Goal Points (The requested single lines) ---
scatter!(p, [START_STATE[1]], [START_STATE[2]], c=:cyan, ms=8, label="Start Point")
scatter!(p, [TARGET_STATE[1]], [TARGET_STATE[2]], c=:red, ms=8, label="Goal Point")
display(p)

# Race Car dynamic with Kinematic Unicycle with Lateral Slip
function dynamics(s, u, n_slip, n_steer, dt, V)
    x, y, th = s

    v_body_x = V
    v_body_y = n_slip 
    
    # 2. Rotate to Global Frame
    # X_dot = vx*cos(th) - vy*sin(th)
    # Y_dot = vx*sin(th) + vy*cos(th)
    dx = v_body_x * cos(th) - v_body_y * sin(th)
    dy = v_body_x * sin(th) + v_body_y * cos(th)
    
    # 3. Heading Update
    # Apply steering control + steering noise
    dth = u + n_steer
    
    # 4. Integrate (Euler Integration)
    xn = x + dx * dt
    yn = y + dy * dt
    thn = th + dth * dt
    
    # 5. Normalize Angle to [0, 2pi)
    thn = mod(thn, 2*pi)
    
    return [xn, yn, thn]
end

# Physical Properties of Race Car
V          = 0.50   # Forward Speed (m/s)
dt         = 0.20   # Time Step (s)
u_limits   = [-pi/2, pi/2] # Steering Limits (rad/s)
nactions = 21
u_actions = collect(LinRange(u_limits[1], u_limits[2], nactions))

# - Bound = 0.60: The drift will never exceed 0.6 m/s (physical limit of tires).
SIGMA_SLIP = 0.3
BOUND_SLIP = 0.4
DIST_SLIP  = Truncated(Normal(0.0, SIGMA_SLIP), -BOUND_SLIP, BOUND_SLIP)

# 2. Steering Error (Control Precision)
# - Sigma = 0.10: Small jitter in the steering wheel.
# - Bound = 0.30: Error is capped at +/- 0.3 rad/s.
SIGMA_STEER = 0.3
BOUND_STEER = 0.35
DIST_STEER  = Truncated(Normal(0.0, SIGMA_STEER), -BOUND_STEER, BOUND_STEER)

############################################################################## Simulation Example to check the car's behaviour 
# x, y, th = 0.0, 0.0, 0.0
# u_cmd    = 2.0  # Constant Left Turn (2.0 rad/s)

# traj_x = [x]
# traj_y = [y]

# Run for 100 steps (10 seconds)
# for t in 1:100
#     global x, y, th
    
#     # 1. Generate Noise
#     n_slip  = rand(DIST_SLIP)
#     n_steer = rand(DIST_STEER)
    
#     # 2. Pack state
#     state_vector = [x, y, th]
    
#     # dynamics(s, u, n_slip, n_steer, dt, V)
#     next_s = dynamics(state_vector, u_cmd, n_slip, n_steer, dt, V)
    
#     # 4. Update x, y, th
#     x, y, th = next_s
    
#     push!(traj_x, x)
#     push!(traj_y, y)
# end

# # === 4. PLOT ===
# p=plot(traj_x, traj_y, 
#     marker=:circle, 
#     markersize=2, 
#     linealpha=0.5,
#     label="Car Path", 
#     title="Stochastic Dynamics Test (Constant Turn)",
#     aspect_ratio=:equal,
#     xlabel="X (m)", ylabel="Y (m)"
# )

# Now try the car's behaviour in the discrete domain 
# creating tree 
# Create 3D States Array [x, y, th]
states_3d = [[x, y, th] for y in y_grid for x in x_grid for th in th_grid]

# Convert vector of vectors to matrix (3 x N) for KDTree
states_matrix = hcat(states_3d...)
tree = KDTree(states_matrix)


# traj_x_grid = [x]
# traj_y_grid = [y]
# # Run for 60 steps
# for t in 1:100
#     global x, y, th
    
#     # 1. Generate Noise
#     n_slip  = rand(DIST_SLIP)
#     n_steer = rand(DIST_STEER)

#     state_vector = [x, y, th]
#     # A. Continuous Update
#     x, y, th = dynamics(state_vector, u_cmd, n_slip, n_steer, dt, V)

#     # B. KNN Mapping (Find Nearest Grid Point)
#     # Query tree for nearest neighbor (k=1)
#     idxs, dists = knn(tree, [x, y, th], 1)
#     nearest_idx = idxs[1]
    
#     # Retrieve Grid State
#     grid_state = states_3d[nearest_idx]
#     gx, gy = grid_state[1], grid_state[2]
    
#     push!(traj_x_grid, gx)
#     push!(traj_y_grid, gy)
# end

# # Overlay Grid Points
# scatter!(p, traj_x_grid, traj_y_grid, 
#     marker=:square, markersize=3, color=:red, label="Nearest Grid Point (KNN)"
# )


# Creating Transition Matrix 
nstates = length(states_3d)
T = [spzeros(nactions, nstates) for _ in 1:nstates]



N_SLIP_SAMPLES  = 10
N_STEER_SAMPLES = 10
TOTAL_SAMPLES   = N_SLIP_SAMPLES * N_STEER_SAMPLES

slip_samples  = rand(DIST_SLIP, N_SLIP_SAMPLES)
steer_samples = rand(DIST_STEER, N_STEER_SAMPLES)


@time begin
    for is in 1:nstates
        s = states_3d[is]
        
        # Optimization: If state is already Unsafe (Wall/Grass), it is Absorbing.
        if !is_safe(s[1], s[2])
            for ia in 1:nactions
                T[is][ia, is] = 1.0 # Stay in wall forever
            end
            continue
        end

        # For valid states, simulate dynamics
        for ia in 1:nactions
            u = u_actions[ia]
            
            # Iterate through both noise loops
            for n_slip in slip_samples
                for n_steer in steer_samples
                    
                    # CALL DYNAMICS FUNCTION
                    # Use the shared function to guarantee physics consistency
                    next_s_cont = dynamics(s, u, n_slip, n_steer, dt, V)

                    # Find nearest discrete state
                    idxs, _ = knn(tree, next_s_cont, 1)
                    next_idx = idxs[1]
                    
                    # Add probability (Uniform weight for all sample combinations)
                    T[is][ia, next_idx] += 1.0 / TOTAL_SAMPLES
                end
            end
        end
        
        # Progress Indicator (Optional)
        if is % 5000 == 0
            print(".")
        end
    end
end


## creating reward function (indicator function)
function reward(s)
    return is_safe(s[1], s[2]) ? 1.0 : 0.0
end
r_vec = [reward(s) for s in states_3d]

#LP for safety calculation
model = Model(Mosek.Optimizer)
set_optimizer_attribute(model, "MSK_IPAR_INTPNT_BASIS", 0)

@variable(model, 0 <= g[1:nstates] <= 1)
@variable(model, h[1:nstates] >= 0)

alpha_dist = ones(nstates) ./ nstates
@objective(model, Min, sum(alpha_dist[s] * g[s] for s in 1:nstates))

@time begin
    for s in 1:nstates
        # OPTIMIZATION:
        # If a state is Unsafe (r=0), it is physically impossible to be safe there.
        # Force g[s] = 0 and skip constraint generation for this state.
        # This significantly reduces the solver load.
        if r_vec[s] == 0.0
            @constraint(model, g[s] == 0.0)
            continue
        end

        # For Safe states, enforce Bellman Inequalities for all actions
        for a in 1:nactions
            # Get transition probabilities P(s' | s, a)
            # T[s] is a sparse matrix (actions x states)
            P_row = T[s][a, :]
            
            # Constraint 1: Gain Stability
            # g(s) >= Sum[ P(s'|s,a) * g(s') ]
            @constraint(model, g[s] >= dot(P_row, g))
            
            # Constraint 2: Bias/Transient Consistency
            # g(s) + h(s) >= r(s,a) + Sum[ P(s'|s,a) * h(s') ]
            # Since we checked r_vec[s] != 0, we know r(s)=1.0 here.
            @constraint(model, g[s] + h[s] >= 1.0 + dot(P_row, h))
        end
    end
end

optimize!(model)
println("Objective value = ", objective_value(model))
println("Constraint Set Proportion: $(round(count(s -> is_safe(s[1], s[2]), states_3d) / nstates * 100, digits=2))%")

status = termination_status(model)
println("Solution found with status: $status")

g_opt = value.(g)
h_opt = value.(h)
max_safety = maximum(g_opt)

println("Max Safety Value: ", round(max_safety))


# Calculating Q function 
Q_safe = zeros(Float64, nstates, nactions)
for s in 1:nstates
    # If the current state is already unsafe (wall), Q is 0 for all actions
    if !is_safe(states_3d[s][1], states_3d[s][2])
        Q_safe[s, :] .= 0.0
        continue
    end
    
    # Calculate expected future safety for each action
    # T[s] is a sparse matrix (actions x next_states)
    Q_safe[s, :] = T[s] * g_opt
end

# === PROPORTION OF SAFE ACTIONS (99%, 70%, 50%) ===
println("Proportion >= 100% Safe: $(round(count(Q_safe .>= 1.0) / length(Q_safe) * 100, digits=4))%")
println("Proportion >= 70% Safe: $(round(count(Q_safe .>= 0.70) / length(Q_safe) * 100, digits=4))%")
println("Proportion >= 50% Safe: $(round(count(Q_safe .>= 0.50) / length(Q_safe) * 100, digits=4))%")


# ploting Q value function 
heatmap(Q_safe.>= 0.75)

# Optional Solution 
# --- 3.2 Stage Cost ---
Q_LQR = Diagonal([1.0, 1.0, 0.0]) 
R_LQR = 1.0 

function get_stage_cost(s_idx, u_idx)
    s = states_3d[s_idx]
    diff = s .- TARGET_STATE
    return (dot(diff, Q_LQR * diff) + u_actions[u_idx]^2 * R_LQR) * dt
end

C_matrix = zeros(Float64, nstates, nactions)
for s in 1:nstates
    for a in 1:nactions
        C_matrix[s, a] = get_stage_cost(s, a)
    end
end

# --- 3.3 Vectorized VI Prep ---
P_big_rows = Int[]; P_big_cols = Int[]; P_big_vals = Float64[]
for s in 1:nstates
    T_s = T[s]
    for a in 1:nactions
        cols_a, vals_a = findnz(T_s[a, :]) 
        global_row = s + (a - 1) * nstates
        append!(P_big_rows, fill(global_row, length(cols_a)))
        append!(P_big_cols, cols_a)
        append!(P_big_vals, vals_a)
    end
end
P_big = sparse(P_big_rows, P_big_cols, P_big_vals, nstates * nactions, nstates)
C_vec_base = vec(C_matrix)

# --- 3.4 Constrained VI Function (Returns Q_matrix AND V_vector) ---
function compute_hard_constraint_policy(alpha, name)
    println("\n   >>> Solving for '$name' (Threshold: $alpha)...")
    
    # FILTER: Mask actions
    valid_mask = Q_safe .>= alpha
    
    percent_valid = count(valid_mask) / (nstates * nactions) * 100
    println("       Allowed Actions: $(round(percent_valid, digits=2))%")

    # INITIALIZATION
    V_val = fill(1e6, nstates)
    idxs, _ = knn(tree, TARGET_STATE, 1)
    V_val[idxs[1]] = 0.0
    
    gamma = 0.99
    
    # APPLY COST MASK (Infinity cost for unsafe actions)
    C_vec_local = copy(C_vec_base)
    C_vec_local[.!vec(valid_mask)] .= 1e6 

    # VALUE ITERATION
    @time for k in 1:5000
        V_prev = V_val
        Q_expected = P_big * V_prev
        Q_vec = C_vec_local .+ gamma .* Q_expected
        
        Q_mat = reshape(Q_vec, (nstates, nactions))
        V_val = minimum(Q_mat, dims=2)[:]
        
        if maximum(abs.(V_val - V_prev)) < 1e-2
            println("       Converged in $k iterations.")
            break
        end
    end
    
    # --- FINAL Q-FUNCTION CALCULATION ---
    # We re-calculate Q one last time using the converged V
    # Q(s,a) = C(s,a) + gamma * Sum(P(s'|s,a) * V(s'))
    Q_final_vec = C_vec_local .+ gamma .* (P_big * V_val)
    Q_final_mat = reshape(Q_final_vec, (nstates, nactions))
    
    # Generate Legal Actions List
    valid_actions_lookup = [findall(valid_mask[s, :]) for s in 1:nstates]
    
    # Fallback for Stuck States
    stuck_count = 0
    for s in 1:nstates
        if isempty(valid_actions_lookup[s])
            stuck_count += 1
            best_unsafe = argmax(Q_safe[s, :]) # Least unsafe
            push!(valid_actions_lookup[s], best_unsafe)
        end
    end
    println("       Stuck States (Fixed with Fallback): $stuck_count")
    
    return Q_final_mat, V_val, valid_actions_lookup
end

# COMPUTE POLICIES
println("value iteration runs now")
Q_risky, V_risky, acts_risky = compute_hard_constraint_policy(0.60, "Risky (Shortcut)")
Q_safe, V_safe, acts_safe    = compute_hard_constraint_policy(0.99, "Safe (Detour)")



## Simulation Section 
function is_target(x, y)
    return (x - TARGET_STATE[1])^2 + (y - TARGET_STATE[2])^2 <= TARGET_RADIUS^2
end

function is_obstacle(x, y)
    return (x - Obs_x)^2 + (y - Obs_y)^2 <= Obs_r^2
end

function run_sim_stochastic(Q_pol, acts_pol, seed)
    cx, cy, cth = START_STATE
    traj_x = [cx]; traj_y = [cy]
    status = :timeout
    
    max_steps = 1200
    for t in 1:max_steps
        if !is_safe(cx, cy); status = :crash; break; end
        if is_target(cx, cy); status = :success; break; end

        idxs, _ = knn(tree, [cx, cy, mod(cth, 2*pi)], 1)
        s_idx = idxs[1]
        
        legal = acts_pol[s_idx]
        if isempty(legal); status = :stuck; break; end
        
        # GREEDY SELECTION USING PRE-CALCULATED Q MATRIX
        # We just pick the legal action with the minimum Q value
        best_u_idx = legal[argmin(Q_pol[s_idx, legal])]
        best_u = u_actions[best_u_idx]
        
        n_slip = rand(DIST_SLIP); n_steer = rand(DIST_STEER)
        nx, ny, nth = dynamics([cx, cy, cth], best_u, n_slip, n_steer, dt, V)
        cx, cy, cth = nx, ny, nth
        push!(traj_x, cx); push!(traj_y, cy)
    end
    return traj_x, traj_y, status
end

NSIMS = 100
TARGET_RADIUS= 0.30 
trajs_risky = [run_sim_stochastic(Q_risky, acts_risky, i) for i in 1:NSIMS]
trajs_safe  = [run_sim_stochastic(Q_safe, acts_safe, i+NSIMS) for i in 1:NSIMS]


# Stats
succ_risky = count(t -> t[3] == :success, trajs_risky)
crash_risky = count(t -> t[3] == :crash, trajs_risky)
stuck_risky = count(t -> t[3] == :stuck, trajs_risky)
time_risky  = count(t -> t[3] == :timeout, trajs_risky)

succ_safe = count(t -> t[3] == :success, trajs_safe)
crash_safe = count(t -> t[3] == :crash, trajs_safe)
stuck_safe = count(t -> t[3] == :stuck, trajs_safe)
time_safe  = count(t -> t[3] == :timeout, trajs_safe)

println("\n[SIMULATION RESULTS - FULL BREAKDOWN]")
println("Risky (0.60):")
println("   Success: $succ_risky")
println("   Crash:   $crash_risky")
println("   Stuck:   $stuck_risky  <-- Likely here")
println("   Timeout: $time_risky")

println("Safe (0.99):")
println("   Success: $succ_safe")
println("   Crash:   $crash_safe")
println("   Stuck:   $stuck_safe   <-- Likely here")
println("   Timeout: $time_safe")

# Plot
heatmap_data = [is_safe(x, y) ? 1.0 : 0.0 for y in y_grid, x in x_grid]
p = heatmap(x_grid, y_grid, heatmap_data, c=:grays, aspect_ratio=:equal, title="Policy Comparison")

for (tx, ty, s) in trajs_risky
    plot!(p, tx, ty, c=:orange, alpha=0.5, label="")
end
for (tx, ty, s) in trajs_safe
    plot!(p, tx, ty, c=:magenta, alpha=0.7, label="")
end

scatter!(p, [START_STATE[1]], [START_STATE[2]], c=:cyan, ms=8, label="Start")
scatter!(p, [TARGET_STATE[1]], [TARGET_STATE[2]], c=:red, ms=8, label="Target")
display(p)

# # 1. Reshape 1D Value Vector into 3D Grid (Theta x X x Y)
# # Dimensions based on your grid: (32, 51, 51)
# V_3d_risky = reshape(V_risky, (num_points_th, num_points_x, num_points_y))
# V_3d_safe  = reshape(V_safe,  (num_points_th, num_points_x, num_points_y))

# # 2. Project to 2D (Take Minimum Cost across all angles)
# # We want the best possible cost for being at (x,y) regardless of heading
# V_2d_risky = dropdims(minimum(V_3d_risky, dims=1), dims=1)
# V_2d_safe  = dropdims(minimum(V_3d_safe,  dims=1), dims=1)

# # 3. Mask Unsafe Areas (Show ONLY Constraint Set)
# # Iterate through the grid and set V=NaN if !is_safe(x,y)
# for (ix, x) in enumerate(x_grid)
#     for (iy, y) in enumerate(y_grid)
#         if !is_safe(x, y)
#             V_2d_risky[ix, iy] = NaN
#             V_2d_safe[ix, iy]  = NaN
#         end
#     end
# end

# # 4. Handle Unreachable States within Safe Set
# # Even inside the track, some states might be unreachable (1e6). Set to NaN or Max for contrast.
# V_plot_risky = replace(V_2d_risky, 1e6 => NaN)
# V_plot_safe  = replace(V_2d_safe,  1e6 => NaN)

# # 5. Plot Heatmaps
# p1 = heatmap(x_grid, y_grid, V_plot_risky', 
#     c=:viridis, title="Risky Policy Value (Cost-to-Go)", xlabel="X", ylabel="Y",
#     xlims=(-3,3), ylims=(-3,3), aspect_ratio=:equal)

# p2 = heatmap(x_grid, y_grid, V_plot_safe', 
#     c=:viridis, title="Safe Policy Value (Cost-to-Go)", xlabel="X", ylabel="Y",
#     xlims=(-3,3), ylims=(-3,3), aspect_ratio=:equal)

# display(plot(p1, p2, layout=(2,1), size=(600, 900)))