# # Safe Optimal simulation for dubins path. 
# # Written by Saber Omidi 

# # Import required libraries
# using Random
# using Distributions
# using NearestNeighbors
# using SparseArrays
# using LinearAlgebra
# using JuMP
# using Mosek
# using MosekTools
# using Plots
# using Dates
# using DelimitedFiles

# # Constant Seed for Randomness 
# Random.seed!(42)

# # Grid Size  a rectangular x,y and theta for orientation
# num_points_x = 55
# num_points_y = 55
# num_points_th = 32

# ## x y range    
# # Rectangular area: 6m long (x) by 4m wide (y)
# x_grid = collect(LinRange(-3.0, 3.0, num_points_x))
# y_grid = collect(LinRange(-3.0, 3.0, num_points_y))
# th_grid = collect(LinRange(0, 2π, num_points_th))

# # Race track as a constrained set.
# function constraint_set(x_grid, y_grid)
#     cx = Float64[]
#     cy = Float64[]
    
#     for x in x_grid
#         for y in y_grid
#             # Geometry Logic: Distance from the central "skeleton" line
#             dx_rect = max(abs(x) - L_straight/2.0, 0.0)
#             dy_rect = abs(y - y_c)
#             d_center = sqrt(dx_rect^2 + dy_rect^2)

#             # If inside the track boundaries, add to the set
#             if (d_center <= R_outer) && (d_center >= R_inner)
#                 push!(cx, x)
#                 push!(cy, y)
#             end
#         end
#     end
    
#     return cx, cy
# end

# L_straight = 2.0 
# R_outer    = 1.80    
# R_inner    = 0.90    
# xc_left    = -L_straight / 2.0
# xc_right   = L_straight / 2.0
# y_c        = 0.0

# ## plot to test the constraint set 
# # 1. Get the set of valid points
# cx, cy = constraint_set(x_grid, y_grid)

# # # 2. Plot
# # p=scatter(cx, cy, 
# #     aspect_ratio=:equal, 
# #     markersize=1.5, 
# #     color=:black,
# #     legend=false, 
# #     title="Constraint Set: Race Track"
# # )

# # Adding obstacles in the race track to make risky descions 
# # === OBSTACLE ===
# Obs_x = 0.0
# Obs_y = 1.1 # Fixed location
# Obs_r = 0.4

# # === OBSTACLE FUNCTION ===
# function obstacle_set(x_grid, y_grid)
#     ox, oy = Float64[], Float64[]
#     for x in x_grid, y in y_grid
#         if (x - Obs_x)^2 + (y - Obs_y)^2 <= Obs_r^2
#             push!(ox, x); push!(oy, y)
#         end
#     end
#     return ox, oy
# end

# ox, oy = obstacle_set(x_grid, y_grid)

# # Ploting to check to Obstacle location 
# #p = scatter(cx, cy, c=:black, ms=1.5, aspect_ratio=:equal, label="Track", title="Track & Obstacle")

# # # 3. Scatter Plot Obstacle (Red)
# # scatter!(p, ox, oy, c=:red, ms=1.5, label="Obstacle")
# # display(p)


# # Is safe function returns bolean one or zero for those are inside the constraint set and they are not inside the obstable. 
# function is_safe(x, y)
#     # 1. Track Boundary Check
#     dx = max(abs(x) - 1.0, 0.0); dy = abs(y)
#     d = sqrt(dx^2 + dy^2)
#     in_track = (d <= R_outer) && (d >= R_inner)
    
#     # 2. Obstacle Check
#     in_obs = (x - Obs_x)^2 + (y - Obs_y)^2 <= Obs_r^2
    
#     return in_track && !in_obs
# end 

# # ploting to check zero and one  bolean 
# heatmap_data = [is_safe(x, y) ? 1.0 : 0.0 for y in y_grid, x in x_grid]

# p = heatmap(x_grid, y_grid, heatmap_data, 
#     c=:grays, 
#     aspect_ratio=:equal, 
#     title="Safety Map: 1=Safe (White), 0=Unsafe (Black)",
#     xlabel="X (m)", ylabel="Y (m)"
# )
# START_STATE  = [-2.5, -0.15, pi/16] 
# TARGET_STATE = [1.5, 1.35, 0] 

# # --- 4.2 Adding Start and Goal Points (The requested single lines) ---
# scatter!(p, [START_STATE[1]], [START_STATE[2]], c=:cyan, ms=8, label="Start Point")
# scatter!(p, [TARGET_STATE[1]], [TARGET_STATE[2]], c=:red, ms=8, label="Goal Point")
# display(p)

# # Race Car dynamic with Kinematic Unicycle with Lateral Slip
# function dynamics(s, u, n_slip, n_steer, dt, V)
#     x, y, th = s

#     v_body_x = V
#     v_body_y = n_slip 
    
#     # 2. Rotate to Global Frame
#     # X_dot = vx*cos(th) - vy*sin(th)
#     # Y_dot = vx*sin(th) + vy*cos(th)
#     dx = v_body_x * cos(th) - v_body_y * sin(th)
#     dy = v_body_x * sin(th) + v_body_y * cos(th)
    
#     # 3. Heading Update
#     # Apply steering control + steering noise
#     dth = u + n_steer
    
#     # 4. Integrate (Euler Integration)
#     xn = x + dx * dt
#     yn = y + dy * dt
#     thn = th + dth * dt
    
#     # 5. Normalize Angle to [0, 2pi)
#     thn = mod(thn, 2*pi)
    
#     return [xn, yn, thn]
# end

# # Physical Properties of Race Car
# V          = 0.50   # Forward Speed (m/s)
# dt         = 0.20   # Time Step (s)
# u_limits   = [-pi/2, pi/2] # Steering Limits (rad/s)
# nactions = 21
# u_actions = collect(LinRange(u_limits[1], u_limits[2], nactions))

# # - Bound = 0.60: The drift will never exceed 0.6 m/s (physical limit of tires).
# SIGMA_SLIP = 0.3
# BOUND_SLIP = 0.4
# DIST_SLIP  = Truncated(Normal(0.0, SIGMA_SLIP), -BOUND_SLIP, BOUND_SLIP)

# # 2. Steering Error (Control Precision)
# # - Sigma = 0.10: Small jitter in the steering wheel.
# # - Bound = 0.30: Error is capped at +/- 0.3 rad/s.
# SIGMA_STEER = 0.3
# BOUND_STEER = 0.35
# DIST_STEER  = Truncated(Normal(0.0, SIGMA_STEER), -BOUND_STEER, BOUND_STEER)

# ############################################################################## Simulation Example to check the car's behaviour 
# # x, y, th = 0.0, 0.0, 0.0
# # u_cmd    = 2.0  # Constant Left Turn (2.0 rad/s)

# # traj_x = [x]
# # traj_y = [y]

# # Run for 100 steps (10 seconds)
# # for t in 1:100
# #     global x, y, th
    
# #     # 1. Generate Noise
# #     n_slip  = rand(DIST_SLIP)
# #     n_steer = rand(DIST_STEER)
    
# #     # 2. Pack state
# #     state_vector = [x, y, th]
    
# #     # dynamics(s, u, n_slip, n_steer, dt, V)
# #     next_s = dynamics(state_vector, u_cmd, n_slip, n_steer, dt, V)
    
# #     # 4. Update x, y, th
# #     x, y, th = next_s
    
# #     push!(traj_x, x)
# #     push!(traj_y, y)
# # end

# # # === 4. PLOT ===
# # p=plot(traj_x, traj_y, 
# #     marker=:circle, 
# #     markersize=2, 
# #     linealpha=0.5,
# #     label="Car Path", 
# #     title="Stochastic Dynamics Test (Constant Turn)",
# #     aspect_ratio=:equal,
# #     xlabel="X (m)", ylabel="Y (m)"
# # )

# # Now try the car's behaviour in the discrete domain 
# # creating tree 
# # Create 3D States Array [x, y, th]
# states_3d = [[x, y, th] for y in y_grid for x in x_grid for th in th_grid]

# # Convert vector of vectors to matrix (3 x N) for KDTree
# states_matrix = hcat(states_3d...)
# tree = KDTree(states_matrix)


# # traj_x_grid = [x]
# # traj_y_grid = [y]
# # # Run for 60 steps
# # for t in 1:100
# #     global x, y, th
    
# #     # 1. Generate Noise
# #     n_slip  = rand(DIST_SLIP)
# #     n_steer = rand(DIST_STEER)

# #     state_vector = [x, y, th]
# #     # A. Continuous Update
# #     x, y, th = dynamics(state_vector, u_cmd, n_slip, n_steer, dt, V)

# #     # B. KNN Mapping (Find Nearest Grid Point)
# #     # Query tree for nearest neighbor (k=1)
# #     idxs, dists = knn(tree, [x, y, th], 1)
# #     nearest_idx = idxs[1]
    
# #     # Retrieve Grid State
# #     grid_state = states_3d[nearest_idx]
# #     gx, gy = grid_state[1], grid_state[2]
    
# #     push!(traj_x_grid, gx)
# #     push!(traj_y_grid, gy)
# # end

# # # Overlay Grid Points
# # scatter!(p, traj_x_grid, traj_y_grid, 
# #     marker=:square, markersize=3, color=:red, label="Nearest Grid Point (KNN)"
# # )


# # Creating Transition Matrix 
# nstates = length(states_3d)
# T = [spzeros(nactions, nstates) for _ in 1:nstates]



# N_SLIP_SAMPLES  = 10
# N_STEER_SAMPLES = 10
# TOTAL_SAMPLES   = N_SLIP_SAMPLES * N_STEER_SAMPLES

# slip_samples  = rand(DIST_SLIP, N_SLIP_SAMPLES)
# steer_samples = rand(DIST_STEER, N_STEER_SAMPLES)


# @time begin
#     for is in 1:nstates
#         s = states_3d[is]
        
#         # Optimization: If state is already Unsafe (Wall/Grass), it is Absorbing.
#         if !is_safe(s[1], s[2])
#             for ia in 1:nactions
#                 T[is][ia, is] = 1.0 # Stay in wall forever
#             end
#             continue
#         end

#         # For valid states, simulate dynamics
#         for ia in 1:nactions
#             u = u_actions[ia]
            
#             # Iterate through both noise loops
#             for n_slip in slip_samples
#                 for n_steer in steer_samples
                    
#                     # CALL DYNAMICS FUNCTION
#                     # Use the shared function to guarantee physics consistency
#                     next_s_cont = dynamics(s, u, n_slip, n_steer, dt, V)

#                     # Find nearest discrete state
#                     idxs, _ = knn(tree, next_s_cont, 1)
#                     next_idx = idxs[1]
                    
#                     # Add probability (Uniform weight for all sample combinations)
#                     T[is][ia, next_idx] += 1.0 / TOTAL_SAMPLES
#                 end
#             end
#         end
        
#         # Progress Indicator (Optional)
#         if is % 5000 == 0
#             print(".")
#         end
#     end
# end


# ## creating reward function (indicator function)
# function reward(s)
#     return is_safe(s[1], s[2]) ? 1.0 : 0.0
# end
# r_vec = [reward(s) for s in states_3d]

# #LP for safety calculation
# model = Model(Mosek.Optimizer)
# set_optimizer_attribute(model, "MSK_IPAR_INTPNT_BASIS", 0)

# @variable(model, 0 <= g[1:nstates] <= 1)
# @variable(model, h[1:nstates] >= 0)

# alpha_dist = ones(nstates) ./ nstates
# @objective(model, Min, sum(alpha_dist[s] * g[s] for s in 1:nstates))

# @time begin
#     for s in 1:nstates
#         # OPTIMIZATION:
#         # If a state is Unsafe (r=0), it is physically impossible to be safe there.
#         # Force g[s] = 0 and skip constraint generation for this state.
#         # This significantly reduces the solver load.
#         if r_vec[s] == 0.0
#             @constraint(model, g[s] == 0.0)
#             continue
#         end

#         # For Safe states, enforce Bellman Inequalities for all actions
#         for a in 1:nactions
#             # Get transition probabilities P(s' | s, a)
#             # T[s] is a sparse matrix (actions x states)
#             P_row = T[s][a, :]
            
#             # Constraint 1: Gain Stability
#             # g(s) >= Sum[ P(s'|s,a) * g(s') ]
#             @constraint(model, g[s] >= dot(P_row, g))
            
#             # Constraint 2: Bias/Transient Consistency
#             # g(s) + h(s) >= r(s,a) + Sum[ P(s'|s,a) * h(s') ]
#             # Since we checked r_vec[s] != 0, we know r(s)=1.0 here.
#             @constraint(model, g[s] + h[s] >= 1.0 + dot(P_row, h))
#         end
#     end
# end

# optimize!(model)
# println("Objective value = ", objective_value(model))
# println("Constraint Set Proportion: $(round(count(s -> is_safe(s[1], s[2]), states_3d) / nstates * 100, digits=2))%")

# status = termination_status(model)
# println("Solution found with status: $status")

# g_opt = value.(g)
# h_opt = value.(h)
# max_safety = maximum(g_opt)

# println("Max Safety Value: ", round(max_safety))


# # Calculating Q function 
# Q_safe = zeros(Float64, nstates, nactions)
# for s in 1:nstates
#     # If the current state is already unsafe (wall), Q is 0 for all actions
#     if !is_safe(states_3d[s][1], states_3d[s][2])
#         Q_safe[s, :] .= 0.0
#         continue
#     end
    
#     # Calculate expected future safety for each action
#     # T[s] is a sparse matrix (actions x next_states)
#     Q_safe[s, :] = T[s] * g_opt
# end

# # === PROPORTION OF SAFE ACTIONS (99%, 70%, 50%) ===
# println("Proportion >= 100% Safe: $(round(count(Q_safe .>= 1.0) / length(Q_safe) * 100, digits=4))%")
# println("Proportion >= 70% Safe: $(round(count(Q_safe .>= 0.70) / length(Q_safe) * 100, digits=4))%")
# println("Proportion >= 50% Safe: $(round(count(Q_safe .>= 0.50) / length(Q_safe) * 100, digits=4))%")


# # ploting Q value function 
# heatmap(Q_safe.>= 0.75)

# # Optional Solution 
# # --- 3.2 Stage Cost ---
# Q_LQR = 10*Diagonal([1.0, 1.0, 0.0]) 
# R_LQR = 0.1 

# function get_stage_cost(s_idx, u_idx)
#     s = states_3d[s_idx]
#     diff = s .- TARGET_STATE
#     return (dot(diff, Q_LQR * diff) + u_actions[u_idx]^2 * R_LQR) * dt
# end

# C_matrix = zeros(Float64, nstates, nactions)
# for s in 1:nstates
#     for a in 1:nactions
#         C_matrix[s, a] = get_stage_cost(s, a)
#     end
# end

# # --- 3.3 Vectorized VI Prep ---
# P_big_rows = Int[]; P_big_cols = Int[]; P_big_vals = Float64[]
# for s in 1:nstates
#     T_s = T[s]
#     for a in 1:nactions
#         cols_a, vals_a = findnz(T_s[a, :]) 
#         global_row = s + (a - 1) * nstates
#         append!(P_big_rows, fill(global_row, length(cols_a)))
#         append!(P_big_cols, cols_a)
#         append!(P_big_vals, vals_a)
#     end
# end
# P_big = sparse(P_big_rows, P_big_cols, P_big_vals, nstates * nactions, nstates)
# C_vec_base = vec(C_matrix)

# # --- 3.4 Constrained VI Function (Returns Q_matrix AND V_vector) ---
# function compute_hard_constraint_policy(alpha, name)
#     println("\n   >>> Solving for '$name' (Threshold: $alpha)...")
    
#     # FILTER: Mask actions
#     valid_mask = Q_safe .>= alpha
    
#     percent_valid = count(valid_mask) / (nstates * nactions) * 100
#     println("       Allowed Actions: $(round(percent_valid, digits=2))%")

#     # INITIALIZATION
#     V_val = fill(1e6, nstates)
#     idxs, _ = knn(tree, TARGET_STATE, 1)
#     V_val[idxs[1]] = 0.0
    
#     gamma = 0.99
    
#     # APPLY COST MASK (Infinity cost for unsafe actions)
#     C_vec_local = copy(C_vec_base)
#     C_vec_local[.!vec(valid_mask)] .= 1e6 

#     # VALUE ITERATION
#     @time for k in 1:5000
#         V_prev = V_val
#         Q_expected = P_big * V_prev
#         Q_vec = C_vec_local .+ gamma .* Q_expected
        
#         Q_mat = reshape(Q_vec, (nstates, nactions))
#         V_val = minimum(Q_mat, dims=2)[:]
        
#         if maximum(abs.(V_val - V_prev)) < 1e-2
#             println("       Converged in $k iterations.")
#             break
#         end
#     end
    
#     # --- FINAL Q-FUNCTION CALCULATION ---
#     # We re-calculate Q one last time using the converged V
#     # Q(s,a) = C(s,a) + gamma * Sum(P(s'|s,a) * V(s'))
#     Q_final_vec = C_vec_local .+ gamma .* (P_big * V_val)
#     Q_final_mat = reshape(Q_final_vec, (nstates, nactions))
    
#     # Generate Legal Actions List
#     valid_actions_lookup = [findall(valid_mask[s, :]) for s in 1:nstates]
    
#     # Fallback for Stuck States
#     stuck_count = 0
#     for s in 1:nstates
#         if isempty(valid_actions_lookup[s])
#             stuck_count += 1
#             best_unsafe = argmax(Q_safe[s, :]) # Least unsafe
#             push!(valid_actions_lookup[s], best_unsafe)
#         end
#     end
#     println("       Stuck States (Fixed with Fallback): $stuck_count")
    
#     return Q_final_mat, V_val, valid_actions_lookup
# end

# # COMPUTE POLICIES
# println("value iteration runs now")
# Q_risky, V_risky, acts_risky = compute_hard_constraint_policy(0.60, "Risky (Shortcut)")
# Q_safe, V_safe, acts_safe    = compute_hard_constraint_policy(0.99, "Safe (Detour)")



# ## Simulation Section 
# function is_target(x, y)
#     return (x - TARGET_STATE[1])^2 + (y - TARGET_STATE[2])^2 <= TARGET_RADIUS^2
# end

# function is_obstacle(x, y)
#     return (x - Obs_x)^2 + (y - Obs_y)^2 <= Obs_r^2
# end

# function run_sim_stochastic(Q_pol, acts_pol, seed)
#     cx, cy, cth = START_STATE
#     traj_x = [cx]; traj_y = [cy]
#     status = :timeout
    
#     max_steps = 1200
#     for t in 1:max_steps
#         if !is_safe(cx, cy); status = :crash; break; end
#         if is_target(cx, cy); status = :success; break; end

#         idxs, _ = knn(tree, [cx, cy, mod(cth, 2*pi)], 1)
#         s_idx = idxs[1]
        
#         legal = acts_pol[s_idx]
#         if isempty(legal); status = :stuck; break; end
        
#         # GREEDY SELECTION USING PRE-CALCULATED Q MATRIX
#         # We just pick the legal action with the minimum Q value
#         best_u_idx = legal[argmin(Q_pol[s_idx, legal])]
#         best_u = u_actions[best_u_idx]
        
#         n_slip = rand(DIST_SLIP); n_steer = rand(DIST_STEER)
#         nx, ny, nth = dynamics([cx, cy, cth], best_u, n_slip, n_steer, dt, V)
#         cx, cy, cth = nx, ny, nth
#         push!(traj_x, cx); push!(traj_y, cy)
#     end
#     return traj_x, traj_y, status
# end

# NSIMS = 100
# TARGET_RADIUS= 0.30 
# trajs_risky = [run_sim_stochastic(Q_risky, acts_risky, i) for i in 1:NSIMS]
# trajs_safe  = [run_sim_stochastic(Q_safe, acts_safe, i+NSIMS) for i in 1:NSIMS]


# # Stats
# succ_risky = count(t -> t[3] == :success, trajs_risky)
# crash_risky = count(t -> t[3] == :crash, trajs_risky)
# stuck_risky = count(t -> t[3] == :stuck, trajs_risky)
# time_risky  = count(t -> t[3] == :timeout, trajs_risky)

# succ_safe = count(t -> t[3] == :success, trajs_safe)
# crash_safe = count(t -> t[3] == :crash, trajs_safe)
# stuck_safe = count(t -> t[3] == :stuck, trajs_safe)
# time_safe  = count(t -> t[3] == :timeout, trajs_safe)

# println("\n[SIMULATION RESULTS - FULL BREAKDOWN]")
# println("Risky (0.60):")
# println("   Success: $succ_risky")
# println("   Crash:   $crash_risky")
# println("   Stuck:   $stuck_risky  <-- Likely here")
# println("   Timeout: $time_risky")

# println("Safe (0.99):")
# println("   Success: $succ_safe")
# println("   Crash:   $crash_safe")
# println("   Stuck:   $stuck_safe   <-- Likely here")
# println("   Timeout: $time_safe")

# # Plot
# heatmap_data = [is_safe(x, y) ? 1.0 : 0.0 for y in y_grid, x in x_grid]
# p = heatmap(x_grid, y_grid, heatmap_data, c=:grays, aspect_ratio=:equal, title="Policy Comparison")

# for (tx, ty, s) in trajs_risky
#     plot!(p, tx, ty, c=:orange, alpha=0.5, label="")
# end
# for (tx, ty, s) in trajs_safe
#     plot!(p, tx, ty, c=:magenta, alpha=0.7, label="")
# end

# scatter!(p, [START_STATE[1]], [START_STATE[2]], c=:cyan, ms=8, label="Start")
# scatter!(p, [TARGET_STATE[1]], [TARGET_STATE[2]], c=:red, ms=8, label="Target")
# display(p)






# # =========================================================
# # SECTION 5: SAVING RESULTS (DATA EXPORT)
# # =========================================================
# println("\n[5/5] Saving Results to 'safe-optimal' folder...")

# # 1. Create Directory
# const RESULT_DIR = "safe-optimal"
# if !isdir(RESULT_DIR)
#     mkdir(RESULT_DIR)
#     println("   Created directory: $RESULT_DIR")
# else
#     println("   Directory exists: $RESULT_DIR")
# end

# # 2. Save Scalar/Vector Data (Gain and Value Functions)
# # We save these as raw text files for easy loading in other tools (MATLAB/Python)
# writedlm(joinpath(RESULT_DIR, "safety_values_g.csv"), g_opt)
# writedlm(joinpath(RESULT_DIR, "value_function_risky.csv"), V_risky)
# writedlm(joinpath(RESULT_DIR, "value_function_safe.csv"), V_safe)

# println("Saved Value Functions (g, V_risky, V_safe).")

# # 3. Save Trajectories to CSV
# # Function to format trajectory data into a matrix
# function export_trajectories(trajs, filename)
#     # Header: RunID, Step, X, Y, Status
#     # We will flatten the list of trajectories into a single big CSV
    
#     data_rows = []
    
#     for (run_id, (tx, ty, status)) in enumerate(trajs)
#         status_str = string(status) # Convert symbol to string
#         for step in 1:length(tx)
#             push!(data_rows, [run_id, step, tx[step], ty[step], status_str])
#         end
#     end
    
#     # Convert to Matrix for DelimitedFiles
#     # Note: Julia's writedlm works best with matrices
#     # We construct a matrix of Any to hold numbers and strings
#     matrix_data = permutedims(hcat(data_rows...))
    
#     header = ["RunID" "Step" "X" "Y" "Status"]
#     final_data = vcat(header, matrix_data)
    
#     writedlm(joinpath(RESULT_DIR, filename), final_data, ',')
# end

# # Export Risky Policy Data
# export_trajectories(trajs_risky, "risky_trajectories.csv")
# println("   Saved Risky Trajectories to CSV.")

# # Export Safe Policy Data
# export_trajectories(trajs_safe, "safe_trajectories.csv")
# println("   Saved Safe Trajectories to CSV.")

# println("\nAll Data Saved Successfully.")


# Safe Optimal simulation for Dubins path (Finish Line + Feasibility Check)
# Rewritten by Gemini based on Saber Omidi's Logic

# Safe Optimal simulation for Dubins path (Finish Line + Data Saving)
# Updated by Gemini

# # Import required libraries
# using Random
# using Distributions
# using NearestNeighbors
# using SparseArrays
# using LinearAlgebra
# using JuMP
# using Mosek
# using MosekTools
# using Plots
# using Dates
# using DelimitedFiles  # Required for saving data

# println("--- Starting Simulation Script (Finish Line Version) ---")

# # Constant Seed for Randomness 
# Random.seed!(42)

# # ==========================================
# # 1. GRID & GEOMETRY SETUP
# # ==========================================
# num_points_x = 55
# num_points_y = 55
# num_points_th = 32

# # Rectangular area: 6m long (x) by 4m wide (y)
# x_grid = collect(LinRange(-3.0, 3.0, num_points_x))
# y_grid = collect(LinRange(-3.0, 3.0, num_points_y))
# th_grid = collect(LinRange(0, 2π, num_points_th))

# # Track Geometry Constants
# L_straight = 2.0 
# R_outer    = 1.80    
# R_inner    = 0.90    
# y_c        = 0.0

# # Obstacle Constants
# Obs_x = 0.0
# Obs_y = 1.1 
# Obs_r = 0.4

# # Start Point
# START_STATE  = [-2.5, -0.15, pi/4] 

# # === FINISH LINE DEFINITION ===
# # Finish line is at x = 1.5 on the bottom curve (y < 0)
# # (Corrected to positive 1.5 so the race covers the full track length)
# FINISH_LINE_X = -1.5 
# FINISH_Y_INNER = -sqrt(R_inner^2 - 0.5^2) # approx -0.748
# FINISH_Y_OUTER = -sqrt(R_outer^2 - 0.5^2) # approx -1.729

# # Target state (center of finish line) for heuristic calculations
# TARGET_STATE = [FINISH_LINE_X, (FINISH_Y_INNER + FINISH_Y_OUTER)/2.0, 0] 

# # Function: Is the point inside the track and outside the obstacle?
# function is_safe(x, y)
#     dx_rect = max(abs(x) - L_straight/2.0, 0.0)
#     dy_rect = abs(y - y_c)
#     d_center = sqrt(dx_rect^2 + dy_rect^2)

#     in_track = (d_center <= R_outer) && (d_center >= R_inner)
#     in_obs = (x - Obs_x)^2 + (y - Obs_y)^2 <= Obs_r^2
    
#     return in_track && !in_obs
# end 

# # Function: Has the car crossed the finish line?
# function is_target(x, y)
#     x_hit = abs(x - FINISH_LINE_X) <= 0.20
#     # Ensure we are on the bottom part of the track (y < 0)
#     y_hit = (y <= FINISH_Y_INNER) && (y >= FINISH_Y_OUTER)
#     return x_hit && y_hit
# end

# # ==========================================
# # 2. DYNAMICS & DISCRETIZATION
# # ==========================================
# function dynamics(s, u, n_slip, n_steer, dt, V)
#     x, y, th = s
#     v_body_x = V
#     v_body_y = n_slip 
#     dx = v_body_x * cos(th) - v_body_y * sin(th)
#     dy = v_body_x * sin(th) + v_body_y * cos(th)
#     dth = u + n_steer
#     xn = x + dx * dt
#     yn = y + dy * dt
#     thn = th + dth * dt
#     thn = mod(thn, 2*pi)
#     return [xn, yn, thn]
# end

# # Physical Properties
# V          = 0.50   
# dt         = 0.20   
# u_limits   = [-pi/3, pi/3] # User updated limits
# nactions   = 21
# u_actions  = collect(LinRange(u_limits[1], u_limits[2], nactions))

# # Noise Distributions (User updated noise)
# SIGMA_SLIP = 0.2
# BOUND_SLIP = 0.3
# DIST_SLIP  = Truncated(Normal(0.0, SIGMA_SLIP), -BOUND_SLIP, BOUND_SLIP)

# SIGMA_STEER = 0.2
# BOUND_STEER = 0.25
# DIST_STEER  = Truncated(Normal(0.0, SIGMA_STEER), -BOUND_STEER, BOUND_STEER)

# # Create 3D States Array & KDTree
# println("Building KDTree and State Space...")
# states_3d = [[x, y, th] for y in y_grid for x in x_grid for th in th_grid]
# states_matrix = hcat(states_3d...)
# tree = KDTree(states_matrix)

# # ==========================================
# # 3. TRANSITION MATRIX (T) CALCULATION
# # ==========================================
# println("Calculating Transition Matrix (this may take a moment)...")
# nstates = length(states_3d)
# T = [spzeros(nactions, nstates) for _ in 1:nstates]

# N_SLIP_SAMPLES  = 10
# N_STEER_SAMPLES = 10
# TOTAL_SAMPLES   = N_SLIP_SAMPLES * N_STEER_SAMPLES
# slip_samples  = rand(DIST_SLIP, N_SLIP_SAMPLES)
# steer_samples = rand(DIST_STEER, N_STEER_SAMPLES)

# @time begin
#     for is in 1:nstates
#         s = states_3d[is]
        
#         # If unsafe, state is absorbing
#         if !is_safe(s[1], s[2])
#             for ia in 1:nactions
#                 T[is][ia, is] = 1.0 
#             end
#             continue
#         end

#         # Calculate Dynamics
#         for ia in 1:nactions
#             u = u_actions[ia]
#             for n_slip in slip_samples
#                 for n_steer in steer_samples
#                     next_s_cont = dynamics(s, u, n_slip, n_steer, dt, V)
#                     idxs, _ = knn(tree, next_s_cont, 1)
#                     next_idx = idxs[1]
#                     T[is][ia, next_idx] += 1.0 / TOTAL_SAMPLES
#                 end
#             end
#         end
#         if is % 10000 == 0; print("."); end
#     end
# end
# println("\nTransition Matrix Complete.")

# # ==========================================
# # 4. SAFETY OPTIMIZATION (LP)
# # ==========================================
# function reward(s)
#     return is_safe(s[1], s[2]) ? 1.0 : 0.0
# end
# r_vec = [reward(s) for s in states_3d]

# println("Setting up Safety LP with Mosek...")
# model = Model(Mosek.Optimizer)
# set_optimizer_attribute(model, "MSK_IPAR_INTPNT_BASIS", 0)
# set_optimizer_attribute(model, "MSK_IPAR_LOG", 0) 

# @variable(model, 0 <= g[1:nstates] <= 1)
# @variable(model, h[1:nstates] >= 0)

# alpha_dist = ones(nstates) ./ nstates
# @objective(model, Min, sum(alpha_dist[s] * g[s] for s in 1:nstates))

# @time begin
#     for s in 1:nstates
#         if r_vec[s] == 0.0
#             @constraint(model, g[s] == 0.0)
#             continue
#         end
#         for a in 1:nactions
#             P_row = T[s][a, :]
#             @constraint(model, g[s] >= dot(P_row, g))
#             @constraint(model, g[s] + h[s] >= 1.0 + dot(P_row, h))
#         end
#     end
# end

# optimize!(model)
# g_opt = value.(g)

# # Calculate Q_safe
# Q_safe = zeros(Float64, nstates, nactions)
# for s in 1:nstates
#     if !is_safe(states_3d[s][1], states_3d[s][2])
#         Q_safe[s, :] .= 0.0
#         continue
#     end
#     Q_safe[s, :] = T[s] * g_opt
# end

# # ==========================================
# # 5. VALUE ITERATION (CONSTRAINED)
# # ==========================================
# Q_LQR = 10*Diagonal([1.0, 1.0, 0.0]) 
# R_LQR = 0.1 

# function get_stage_cost(s_idx, u_idx)
#     s = states_3d[s_idx]
#     diff = s .- TARGET_STATE
#     return (dot(diff, Q_LQR * diff) + u_actions[u_idx]^2 * R_LQR) * dt
# end

# C_matrix = zeros(Float64, nstates, nactions)
# for s in 1:nstates
#     for a in 1:nactions
#         C_matrix[s, a] = get_stage_cost(s, a)
#     end
# end

# # Prepare Big Sparse Matrix for VI
# P_big_rows = Int[]; P_big_cols = Int[]; P_big_vals = Float64[]
# for s in 1:nstates
#     T_s = T[s]
#     for a in 1:nactions
#         cols_a, vals_a = findnz(T_s[a, :]) 
#         global_row = s + (a - 1) * nstates
#         append!(P_big_rows, fill(global_row, length(cols_a)))
#         append!(P_big_cols, cols_a)
#         append!(P_big_vals, vals_a)
#     end
# end
# P_big = sparse(P_big_rows, P_big_cols, P_big_vals, nstates * nactions, nstates)
# C_vec_base = vec(C_matrix)

# function compute_hard_constraint_policy(alpha, name)
#     println("\n   >>> Solving for '$name' (Threshold: $alpha)...")
#     valid_mask = Q_safe .>= alpha
    
#     # Initialize V
#     V_val = fill(1e6, nstates)
#     idxs, _ = knn(tree, TARGET_STATE, 1)
#     V_val[idxs[1]] = 0.0
#     gamma = 0.99
    
#     # Apply Mask (Infinite cost for unsafe actions)
#     C_vec_local = copy(C_vec_base)
#     C_vec_local[.!vec(valid_mask)] .= 1e6 

#     # Value Iteration Loop
#     for k in 1:5000
#         V_prev = V_val
#         Q_expected = P_big * V_prev
#         Q_vec = C_vec_local .+ gamma .* Q_expected
#         Q_mat = reshape(Q_vec, (nstates, nactions))
#         V_val = minimum(Q_mat, dims=2)[:]
#         if maximum(abs.(V_val - V_prev)) < 1e-2
#             break
#         end
#     end
    
#     Q_final_vec = C_vec_local .+ gamma .* (P_big * V_val)
#     Q_final_mat = reshape(Q_final_vec, (nstates, nactions))
    
#     valid_actions_lookup = [findall(valid_mask[s, :]) for s in 1:nstates]
#     return Q_final_mat, V_val, valid_actions_lookup
# end

# Q_risky, V_risky, acts_risky = compute_hard_constraint_policy(0.60, "Risky (Shortcut)")
# Q_safe, V_safe, acts_safe    = compute_hard_constraint_policy(0.95, "Safe (Detour)")

# # ==========================================
# # 6. DISCRETE SIMULATION
# # ==========================================
# function sample_next_state(sparse_row)
#     cols, vals = findnz(sparse_row)
#     if isempty(vals); return nothing; end
#     r = rand()
#     cum_sum = 0.0
#     for i in 1:length(vals)
#         cum_sum += vals[i]
#         if r <= cum_sum; return cols[i]; end
#     end
#     return cols[end]
# end

# function run_sim_discrete(Q_pol, acts_pol, seed)
#     Random.seed!(seed)
#     idxs, _ = knn(tree, START_STATE, 1)
#     s_idx = idxs[1]
    
#     cx, cy = states_3d[s_idx][1], states_3d[s_idx][2]
#     traj_x = [cx]; traj_y = [cy]
#     status = :timeout
    
#     for t in 1:1200
#         curr_state = states_3d[s_idx]
#         cx, cy = curr_state[1], curr_state[2]
        
#         if !is_safe(cx, cy); status = :crash; break; end
#         if is_target(cx, cy); status = :success; break; end
        
#         legal = acts_pol[s_idx]
#         if isempty(legal); status = :stuck; break; end
        
#         best_u_idx = legal[argmin(Q_pol[s_idx, legal])]
#         prob_row = T[s_idx][best_u_idx, :]
#         next_idx = sample_next_state(prob_row)
        
#         if isnothing(next_idx); status = :stuck; break; end
#         s_idx = next_idx
#         push!(traj_x, states_3d[s_idx][1])
#         push!(traj_y, states_3d[s_idx][2])
#     end
#     return traj_x, traj_y, status
# end

# println("Starting Simulations...")
# NSIMS = 100
# trajs_risky = [run_sim_discrete(Q_risky, acts_risky, i) for i in 1:NSIMS]
# trajs_safe  = [run_sim_discrete(Q_safe, acts_safe, i+NSIMS) for i in 1:NSIMS]

# succ_risky  = count(t -> t[3] == :success, trajs_risky)
# stuck_risky = count(t -> t[3] == :stuck, trajs_risky)
# crash_risky = count(t -> t[3] == :crash, trajs_risky)

# succ_safe   = count(t -> t[3] == :success, trajs_safe)
# stuck_safe  = count(t -> t[3] == :stuck, trajs_safe)
# crash_safe  = count(t -> t[3] == :crash, trajs_safe)

# println("\n[RESULTS]")
# println("Risky (0.60) | Success: $succ_risky | No Feasible: $stuck_risky | Crash: $crash_risky")
# println("Safe  (0.95) | Success: $succ_safe | No Feasible: $stuck_safe  | Crash: $crash_safe")

# # ==========================================
# # 7. PLOTTING WITH FINISH LINE
# # ==========================================
# println("\nGenerating Plot...")

# # Background
# heatmap_data = [is_safe(x, y) ? 1.0 : 0.0 for y in y_grid, x in x_grid]
# p = heatmap(x_grid, y_grid, heatmap_data, 
#     c=:grays, aspect_ratio=:equal, title="Simulation with Finish Line", legend=false)

# # Trajectories
# for (tx, ty, s) in trajs_risky; plot!(p, tx, ty, c=:orange, alpha=0.5); end
# for (tx, ty, s) in trajs_safe; plot!(p, tx, ty, c=:magenta, alpha=0.5); end

# # Start Point
# scatter!(p, [START_STATE[1]], [START_STATE[2]], c=:cyan, ms=6, label="Start")

# # Draw Finish Line
# plot!(p, 
#     [FINISH_LINE_X, FINISH_LINE_X],   # X coords
#     [FINISH_Y_OUTER, FINISH_Y_INNER], # Y coords
#     color=:red, linewidth=4, label="Finish Line"
# )

# display(p)

# # =========================================================
# # SECTION 8: SAVING RESULTS (DATA EXPORT)
# # =========================================================
# println("\n[5/5] Saving Results to 'DC/safe-optimal' folder...")

# # 1. DEFINE PATH
# # Use @__DIR__ to get the folder where this script is located.
# # Add "DC/safe-optimal" to that path.
# SCRIPT_DIR = @__DIR__
# const RESULT_DIR = joinpath(SCRIPT_DIR, "safe-optimal")

# # mkpath recursively creates directories (e.g., creates DC first, then safe-optimal)
# if !isdir(RESULT_DIR)
#     mkpath(RESULT_DIR)
#     println("   Created directory: $RESULT_DIR")
# else
#     println("   Directory exists: $RESULT_DIR")
# end

# # 2. Save Scalar/Vector Data (Gain and Value Functions)
# writedlm(joinpath(RESULT_DIR, "safety_values_g.csv"), g_opt)
# writedlm(joinpath(RESULT_DIR, "value_function_risky.csv"), V_risky)
# writedlm(joinpath(RESULT_DIR, "value_function_safe.csv"), V_safe)

# println("Saved Value Functions (g, V_risky, V_safe).")

# # 3. Save Trajectories to CSV
# function export_trajectories(trajs, filename)
#     data_rows = []
    
#     for (run_id, (tx, ty, status)) in enumerate(trajs)
#         status_str = string(status) # Convert symbol to string
#         for step in 1:length(tx)
#             push!(data_rows, [run_id, step, tx[step], ty[step], status_str])
#         end
#     end
    
#     if isempty(data_rows)
#         println("   Warning: No trajectory data to save for $filename")
#         return
#     end

#     matrix_data = permutedims(hcat(data_rows...))
#     header = ["RunID" "Step" "X" "Y" "Status"]
#     final_data = vcat(header, matrix_data)
    
#     writedlm(joinpath(RESULT_DIR, filename), final_data, ',')
# end

# # Export Risky Policy Data
# export_trajectories(trajs_risky, "risky_trajectories.csv")
# println("   Saved Risky Trajectories to CSV.")

# # Export Safe Policy Data
# export_trajectories(trajs_safe, "safe_trajectories.csv")
# println("   Saved Safe Trajectories to CSV.")

# println("\nAll Data Saved Successfully.")


# function get_optimal_action_sets(Q_matrix; tol=1e-2)
#     n_states, n_actions = size(Q_matrix)
#     optimal_policy_sets = Vector{Vector{Int}}(undef, n_states)
    
#     for s in 1:n_states
#         # 1. Find the best value (even if it is huge/unsafe)
#         v_star = minimum(Q_matrix[s, :])
        
#         # 2. If the state is TRULY impossible (like essentially infinite everywhere),
#         # we can check if v_star is practically Inf, but usually standard 
#         # masking sets it to 1e6. We want to keep moving even then.
#         # We REMOVED the check: "if v_star > 1e5..."
        
#         optimal_actions = Int[]
#         for a in 1:n_actions
#             # 3. Select all actions that are close to the best value
#             if abs(Q_matrix[s, a] - v_star) < tol
#                 push!(optimal_actions, a)
#             end
#         end
        
#         # Safety Fallback: 
#         # If tolerance failed (rare numerical edge case), force pick argmin
#         if isempty(optimal_actions)
#             push!(optimal_actions, argmin(Q_matrix[s, :]))
#         end
        
#         optimal_policy_sets[s] = optimal_actions
#     end
    
#     return optimal_policy_sets
# end

# # --- Re-run the extraction ---
# println("Extracting Optimal Sets (Correction Applied)...")
# opt_sets_risky = get_optimal_action_sets(Q_risky, tol=1e-2)
# opt_sets_safe  = get_optimal_action_sets(Q_safe,  tol=1e-2)

# # --- Verify Start State ---
# idxs, _ = knn(tree, START_STATE, 1)
# s_idx = idxs[1]
# println("Start State Actions (Risky): ", opt_sets_risky[s_idx])
# println("Start State Value (Risky):   ", minimum(Q_risky[s_idx, :])) 

# # --- Verify Start State for SAFE Policy ---
# idxs, _ = knn(tree, START_STATE, 1)
# s_idx = idxs[1]

# println("Start State Actions (Safe): ", opt_sets_safe[s_idx])
# println("Start State Value (Safe):   ", minimum(Q_safe[s_idx, :]))

# # If Value is > 1e5, it means the start is technically 'unsafe', but now we have actions!
# # =========================================================
# # SECTION 6: OPTIMIZED ANALYTICAL VERIFICATION
# # =========================================================
# println("\n[6/6] Verifying Policy Safety (Optimized Matrix Method)...")

# function compute_policy_safety_gain_optimized(optimal_policy_sets, name)
#     println("   1. Building Global Transition Matrix for '$name'...")
    
#     # We will build a single sparse matrix P_pi (nstates x nstates)
#     # P_pi[s, s'] = Probability of moving s -> s' under the policy
    
#     I_idx = Int[]
#     J_idx = Int[]
#     V_val = Float64[]
    
#     # Pre-allocate estimation (approx 100 non-zeros per state)
#     sizehint!(I_idx, nstates * 100)
#     sizehint!(J_idx, nstates * 100)
#     sizehint!(V_val, nstates * 100)
    
#     # Construct the matrix ONCE
#     # This avoids row-slicing inside the loop
#     for s in 1:nstates
#         best_actions = optimal_policy_sets[s]
        
#         # If no actions or unsafe, we just don't add transitions (row s will be 0)
#         if isempty(best_actions)
#             continue
#         end
        
#         # We need to average the transitions of all optimal actions.
#         # Weight for each action = 1.0 / number of choices
#         weight = 1.0 / length(best_actions)
        
#         # Extract non-zeros from T[s] directly (Fastest method)
#         # T[s] is small (100 entries), so findnz is very fast
#         rows, cols, vals = findnz(T[s])
        
#         for k in 1:length(rows)
#             action_idx = rows[k]
            
#             # Only include transitions for actions that are in our Optimal Set
#             if action_idx in best_actions
#                 push!(I_idx, s)           # Row s (Current State)
#                 push!(J_idx, cols[k])     # Col s' (Next State)
#                 push!(V_val, vals[k] * weight)
#             end
#         end
#     end
    
#     # Create the global sparse matrix
#     P_pi = sparse(I_idx, J_idx, V_val, nstates, nstates)
#     println("      Matrix built with $(nnz(P_pi)) transitions.")

#     # ---------------------------------------------------------
#     # 2. Power Iteration (Now extremely fast matrix-vector mult)
#     # ---------------------------------------------------------
#     println("   2. Solving g = P_pi * g...")
    
#     # Initialize Gain (1.0 everywhere initially)
#     g_val = ones(Float64, nstates)
    
#     # Set strictly unsafe states (walls) to 0.0
#     for s in 1:nstates
#         if !is_safe(states_3d[s][1], states_3d[s][2])
#             g_val[s] = 0.0
#         end
#     end
    
#     max_iter = 3000
#     tol = 1e-6
    
#     @time for k in 1:max_iter
#         g_prev = g_val
        
#         # FAST STEP: One big matrix multiplication
#         # P_pi * g calculates the expected safety for next step
#         g_val = P_pi * g_prev
        
#         # Re-enforce boundary conditions (Walls must stay 0)
#         # (Though naturally P_pi rows for walls should be empty/zero anyway)
#         # Ideally, P_pi handles this, but we can clamp if needed.
        
#         # Convergence Check
#         diff = maximum(abs.(g_val - g_prev))
#         if diff < tol
#             println("      Converged in $k iterations.")
#             break
#         end
#     end
    
#     return g_val
# end

# # --- Run the Optimized Verification ---
# g_verify_risky = compute_policy_safety_gain_optimized(opt_sets_risky, "Risky Policy")
# g_verify_safe  = compute_policy_safety_gain_optimized(opt_sets_safe,  "Safe Policy")

# # Save results
# writedlm(joinpath(RESULT_DIR, "verified_safety_risky.csv"), g_verify_risky)
# writedlm(joinpath(RESULT_DIR, "verified_safety_safe.csv"), g_verify_safe)
# println("Verified safety maps saved.")


# # =========================================================
# # SECTION 9: MONTE CARLO DISTRIBUTION - SAFE POLICY
# # =========================================================
# println("\n[9/10] Generating State Distribution (Safe Policy)...")

# # 1. Setup Simulation Parameters
# mc_policy_Q_safe    = Q_safe      # We are analyzing the SAFE policy
# mc_policy_acts_safe = acts_safe
# N_MC_TRIALS         = 100000        # Number of simulations to run
# MC_HORIZON          = 150         # Max steps per simulation

# # 2. Data Containers
# dist_x_safe = Float64[]
# dist_y_safe = Float64[]

# # Run Simulations
# for i in 1:N_MC_TRIALS
#     # Start with tiny noise for natural spread
#     c_state = START_STATE .+ [rand(Normal(0,0.01)), rand(Normal(0,0.01)), rand(Normal(0,0.01))]
    
#     for t in 1:MC_HORIZON
#         cx, cy, cth = c_state
        
#         # Check Terminal Conditions
#         if !is_safe(cx, cy) || is_target(cx, cy)
#             push!(dist_x_safe, cx); push!(dist_y_safe, cy)
#             break
#         end
        
#         # Store State
#         push!(dist_x_safe, cx); push!(dist_y_safe, cy)

#         # Get Action from Policy
#         idxs, _ = knn(tree, [cx, cy, mod(cth, 2*pi)], 1)
#         s_idx = idxs[1]
#         legal = mc_policy_acts_safe[s_idx]
#         if isempty(legal); break; end # Stuck
        
#         # Greedy choice
#         best_u_idx = legal[argmin(mc_policy_Q_safe[s_idx, legal])]
#         u_cmd = u_actions[best_u_idx]
        
#         # Evolve with Noise
#         n_slip  = rand(DIST_SLIP); n_steer = rand(DIST_STEER)
#         c_state = dynamics(c_state, u_cmd, n_slip, n_steer, dt, V)
#     end
# end

# # 3. Visualization (Safe Policy)
# # Background Track Map (no colorbar for this layer)
# heatmap_data = [is_safe(x, y) ? 1.0 : 0.0 for y in y_grid, x in x_grid]
# p_safe = heatmap(x_grid, y_grid, heatmap_data, 
#     c=:grays, 
#     colorbar=false,
#     legend=false, 
#     title="State Distribution (Safe Policy)\n(Color = Probability Density)",
#     aspect_ratio=:equal,
#     alpha=0.3,
#     xlabel="X (m)", ylabel="Y (m)"
# )

# # Overlay 2D Histogram (Heatmap with Colorbar)
# histogram2d!(p_safe, dist_x_safe, dist_y_safe,
#     bins = (60, 60),       # Higher resolution bins
#     c = :inferno,          # 'inferno' colormap (Black->Purple->Magenta->Yellow)
#     alpha = 0.8,           
#     normalize = :pdf,      # Normalize counts to represent Probability Density
#     colorbar = true        # **Add the Color Bar**
# )

# # Add context markers
# scatter!(p_safe, [START_STATE[1]], [START_STATE[2]], c=:cyan, ms=6, label="Start")
# plot!(p_safe, [FINISH_LINE_X, FINISH_LINE_X], [FINISH_Y_OUTER, FINISH_Y_INNER], 
#     color=:green, linewidth=3, label="Finish")


# # =========================================================
# # SECTION 10: MONTE CARLO DISTRIBUTION - RISKY POLICY
# # =========================================================
# println("\n[10/10] Generating State Distribution (Risky Policy)...")

# # 1. Setup (Switch to Risky Policy)
# mc_policy_Q_risky    = Q_risky
# mc_policy_acts_risky = acts_risky
# # reuse N_MC_TRIALS, MC_HORIZON

# # 2. Data Containers
# dist_x_risky = Float64[]
# dist_y_risky = Float64[]

# # Run Simulations
# for i in 1:N_MC_TRIALS
#     c_state = START_STATE .+ [rand(Normal(0,0.01)), rand(Normal(0,0.01)), rand(Normal(0,0.01))]
    
#     for t in 1:MC_HORIZON
#         cx, cy, cth = c_state
        
#         if !is_safe(cx, cy) || is_target(cx, cy)
#             push!(dist_x_risky, cx); push!(dist_y_risky, cy)
#             break
#         end
#         push!(dist_x_risky, cx); push!(dist_y_risky, cy)

#         # Get Action from RISKY Policy
#         idxs, _ = knn(tree, [cx, cy, mod(cth, 2*pi)], 1)
#         s_idx = idxs[1]
#         legal = mc_policy_acts_risky[s_idx]
#         if isempty(legal); break; end
        
#         best_u_idx = legal[argmin(mc_policy_Q_risky[s_idx, legal])]
#         u_cmd = u_actions[best_u_idx]
        
#         n_slip  = rand(DIST_SLIP); n_steer = rand(DIST_STEER)
#         c_state = dynamics(c_state, u_cmd, n_slip, n_steer, dt, V)
#     end
# end

# # 3. Visualization (Risky Policy)
# # Background Track Map
# p_risky = heatmap(x_grid, y_grid, heatmap_data, 
#     c=:grays, 
#     colorbar=false,
#     legend=false, 
#     title="State Distribution (Risky Policy)\n(Note the tighter turns)",
#     aspect_ratio=:equal,
#     alpha=0.3,
#     xlabel="X (m)" # No ylabel for the second plot to save space
# )

# # Overlay 2D Histogram
# histogram2d!(p_risky, dist_x_risky, dist_y_risky,
#     bins = (60, 60),
#     c = :inferno,
#     alpha = 0.8,           
#     normalize = :pdf,
#     colorbar = true
# )

# # Add context markers
# scatter!(p_risky, [START_STATE[1]], [START_STATE[2]], c=:cyan, ms=6, label="Start")
# plot!(p_risky, [FINISH_LINE_X, FINISH_LINE_X], [FINISH_Y_OUTER, FINISH_Y_INNER], 
#     color=:green, linewidth=3, label="Finish")

# # =========================================================
# # FINAL PLOT DISPLAY (Combine Side-by-Side)
# # =========================================================
# println("Displaying Comparison Plot...")
# # Combine plots in a 1x2 layout and set overall size
# final_plot = plot(p_safe, p_risky, layout=(1, 2), size=(1000, 500))
# display(final_plot)
# println("Analysis Complete.")




# # =========================================================
# # SECTION 11: EXACT DISTRIBUTION VIA LINEAR PROGRAM (DUAL)
# # =========================================================
# println("\n[11/11] Solving Analytical State Distribution (Linear Program)...")

# function solve_occupation_measure_lp(policy_acts, policy_Q, policy_name)
#     println("   >>> Solving Flow LP for: $policy_name")

#     # 1. Build the Induced Transition Matrix (P_pi)
#     # P_pi[i, j] = Probability of moving from state i -> j under policy pi
#     # We need the Transpose P_pi' for the flow constraint: Flow In = Flow Out
    
#     # We build P_pi_T directly (j -> i) for efficiency
#     I_idx = Int[]
#     J_idx = Int[]
#     V_val = Float64[]
    
#     # Pre-allocate
#     sizehint!(I_idx, nstates * 10)
#     sizehint!(J_idx, nstates * 10)
#     sizehint!(V_val, nstates * 10)

#     for s in 1:nstates
#         # Get Action from Policy
#         legal = policy_acts[s]
#         if isempty(legal)
#             # If stuck, no transition out (absorbing)
#             continue
#         end
        
#         # Greedy Policy Selection
#         # (If multiple optimal actions exist, we pick the first one for the flow calc)
#         best_u_idx = legal[argmin(policy_Q[s, legal])]
        
#         # Get transition probabilities for this (state, action) pair
#         # T[s] is [Actions x NextStates]
#         # We want row = best_u_idx
#         rows, cols, vals = findnz(T[s])
        
#         for k in 1:length(rows)
#             if rows[k] == best_u_idx
#                 next_state = cols[k]
#                 prob = vals[k]
                
#                 # Add to Sparse Matrix Construction
#                 # We are building P_transpose: P_ji = Prob(j | i)
#                 # So row = next_state, col = s
#                 push!(I_idx, next_state)
#                 push!(J_idx, s)
#                 push!(V_val, prob)
#             end
#         end
#     end
    
#     # Construct Sparse Matrix P^T
#     P_transpose = sparse(I_idx, J_idx, V_val, nstates, nstates)
    
#     # 2. Setup Linear Program (Flow Balance)
#     # mu = nu + P^T * mu
#     # Rearranged: (I - P^T) * mu = nu
    
#     # Define Source Vector (nu)
#     # The system starts at START_STATE with probability 1.0
#     idxs, _ = knn(tree, START_STATE, 1)
#     start_idx = idxs[1]
    
#     nu = zeros(nstates)
#     nu[start_idx] = 1.0
    
#     # 3. Solve Linear System
#     # Since this is a simple linear equality, we don't technically need an optimizer
#     # like Mosek for the solve, we can use Julia's linear algebra `\` operator.
#     # However, (I - P^T) can be singular if there are perfect cycles.
#     # We add a tiny regularization term (discount factor) for stability if needed.
#     # (I - gamma * P^T) * mu = nu
    
#     gamma_flow = 1 # Discount factor to ensure numerical stability (evaporation)
#     A = I - gamma_flow * P_transpose
    
#     println("      Solving linear system A * mu = nu ...")
#     mu = A \ nu
    
#     # Clean up small negatives from numerical noise
#     mu = max.(mu, 0.0)
    
#     return mu
# end

# # --- RUN FOR BOTH POLICIES ---
# mu_safe = solve_occupation_measure_lp(acts_safe, Q_safe, "Safe Policy")
# mu_risky = solve_occupation_measure_lp(acts_risky, Q_risky, "Risky Policy")

# # =========================================================
# # VISUALIZATION (ANALYTICAL HEATMAP)
# # =========================================================
# println("   Generating Analytical Heatmaps...")

# # Helper to convert vector mu back to 2D grid
# function mu_to_heatmap(mu_vec)
#     grid_map = zeros(num_points_y, num_points_x) # y, x
#     # Sum occupation measures over all orientations for each x,y
#     # states_3d is ordered: y loops outer, x loops middle, th loops inner
#     # Actually, let's just loop to be safe
#     for s in 1:nstates
#         # Determine grid indices from state value (approximate back to grid)
#         # Or better: we know the structure of states_3d construction
#         # But even simpler: let's map mu values to x,y coordinates and let Plots bin them?
#         # No, we want the grid directly.
        
#         # Reverse engineer the index or just use the state coordinate
#         sx, sy = states_3d[s][1], states_3d[s][2]
        
#         # Find x index
#         x_idx = searchsortedfirst(x_grid, sx)
#         y_idx = searchsortedfirst(y_grid, sy)
        
#         if 1 <= x_idx <= num_points_x && 1 <= y_idx <= num_points_y
#             grid_map[y_idx, x_idx] += mu_vec[s]
#         end
#     end
#     return grid_map
# end

# # Convert
# grid_mu_safe = mu_to_heatmap(mu_safe)
# grid_mu_risky = mu_to_heatmap(mu_risky)

# # Plot Safe
# p_lp_safe = heatmap(x_grid, y_grid, grid_mu_safe, 
#     c=:inferno, 
#     title="Analytical Distribution (Safe LP)", 
#     aspect_ratio=:equal,
#     clims=(0, maximum(grid_mu_safe)*0.6) # Saturate slightly to see path
# )
# scatter!(p_lp_safe, [START_STATE[1]], [START_STATE[2]], c=:cyan, label=false)

# # Plot Risky
# p_lp_risky = heatmap(x_grid, y_grid, grid_mu_risky, 
#     c=:inferno, 
#     title="Analytical Distribution (Risky LP)", 
#     aspect_ratio=:equal,
#     clims=(0, maximum(grid_mu_risky)*0.6)
# )
# scatter!(p_lp_risky, [START_STATE[1]], [START_STATE[2]], c=:cyan, label=false)

# # Display Side-by-Side
# final_lp_plot = plot(p_lp_safe, p_lp_risky, layout=(1,2), size=(1000, 450))
# display(final_lp_plot)



# # =========================================================
# # SECTION 10: MONTE CARLO DISTRIBUTION - RISKY POLICY (ICY ROAD)
# # =========================================================
# println("\n[10/10] Generating State Distribution (Risky Policy - High Noise)...")

# # 1. Setup (Risky Policy)
# mc_policy_Q_risky    = Q_risky
# mc_policy_acts_risky = acts_risky
# # Reuse N_MC_TRIALS, MC_HORIZON from Section 9

# # 2. Data Containers
# dist_x_risky = Float64[]
# dist_y_risky = Float64[]

# # === NEW: SIMULATION NOISE MULTIPLIER ===
# # The controller was trained for noise x1.0. 
# # We simulate with noise x3.0 to force collisions for the demo.
# NOISE_FACTOR = 2.0 

# # Run Simulations
# for i in 1:N_MC_TRIALS
#     # Initial state jitter
#     c_state = START_STATE .+ [rand(Normal(0,0.01)), rand(Normal(0,0.01)), rand(Normal(0,0.01))]
    
#     for t in 1:MC_HORIZON
#         cx, cy, cth = c_state
        
#         # Check Terminal Conditions
#         if !is_safe(cx, cy) || is_target(cx, cy)
#             push!(dist_x_risky, cx); push!(dist_y_risky, cy)
#             break
#         end
#         push!(dist_x_risky, cx); push!(dist_y_risky, cy)

#         # Get Action
#         idxs, _ = knn(tree, [cx, cy, mod(cth, 2*pi)], 1)
#         s_idx = idxs[1]
#         legal = mc_policy_acts_risky[s_idx]
#         if isempty(legal); break; end
        
#         # Greedy choice
#         best_u_idx = legal[argmin(mc_policy_Q_risky[s_idx, legal])]
#         u_cmd = u_actions[best_u_idx]
        
#         # === MODIFY PHYSICS HERE ===
#         # Multiply the random noise by the factor
#         n_slip  = rand(DIST_SLIP) * NOISE_FACTOR
#         n_steer = rand(DIST_STEER) * NOISE_FACTOR
        
#         c_state = dynamics(c_state, u_cmd, n_slip, n_steer, dt, V)
#     end
# end

# # 3. Visualization (Risky Policy - Icy)
# # Background Track Map
# p_risky = heatmap(x_grid, y_grid, heatmap_data, 
#     c=:grays, 
#     colorbar=false,
#     legend=false, 
#     title="Risky Policy Distribution\n(Simulated on 3x Noise Surface)",
#     aspect_ratio=:equal,
#     alpha=0.3,
#     xlabel="X (m)"
# )

# # Overlay 2D Histogram
# histogram2d!(p_risky, dist_x_risky, dist_y_risky,
#     bins = (60, 60),
#     c = :inferno,
#     alpha = 0.8,           
#     normalize = :pdf,
#     colorbar = true
# )

# # Context
# scatter!(p_risky, [START_STATE[1]], [START_STATE[2]], c=:cyan, ms=6, label="Start")
# plot!(p_risky, [FINISH_LINE_X, FINISH_LINE_X], [FINISH_Y_OUTER, FINISH_Y_INNER], 
#     color=:green, linewidth=3, label="Finish")

# # Display Comparison
# println("Displaying Comparison Plot...")
# final_plot = plot(p_safe, p_risky, layout=(1, 2), size=(1000, 500))
# display(final_plot)
# println("Analysis Complete.")


## Written By Saber Omidi ##
##########################
ENV["GKSwstype"] = "100"
println("Code Started")
### libraries 
using Random
using Distributions
using NearestNeighbors
using SparseArrays
using LinearAlgebra
using JuMP, Dualization
using Mosek
using MosekTools
using Plots
using Dates
using DelimitedFiles
using Printf



Random.seed!(42) # Constant Seed for Randomness 
plot_test = false    # Swtich for test plotting 
transition_test = false # transition matrix test 

### Geometry and Physics of the problem
# Grid Size  a rectangular x,y and theta for orientation
num_points_x = 61
num_points_y = 61
num_points_th = 40
  
# Rectangular area: 6m long (x) by 4m wide (y)
x_grid = collect(LinRange(-3.0, 3.0, num_points_x))
y_grid = collect(LinRange(-3.0, 3.0, num_points_y))
th_grid = collect(LinRange(0, 2π, num_points_th))


# Race track as a constrained set 
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

# Parameters for constraint set geometry
L_straight = 2.0 
R_outer    = 1.80    
R_inner    = 0.90    
xc_left    = -L_straight / 2.0
xc_right   = L_straight / 2.0
y_c        = 0.0

# plot to test the constraint set 
cx, cy = constraint_set(x_grid, y_grid)

if plot_test 
p=scatter(cx, cy, 
    aspect_ratio=:equal, 
    markersize=1.5, 
    color=:black,
    legend=false, 
    title="Constraint Set: Race Track"
)
end 

# Adding obstacles
Obs_x = 0.0
Obs_y = 1.1 
Obs_r = 0.4

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

if plot_test
p = scatter(cx, cy, c=:black, ms=2.5, aspect_ratio=:equal, label="Track", title="Track & Obstacle")
scatter!(p, ox, oy, c=:red, ms=3.5, label="Obstacle")
end 


# Is safe function returns bolean one or zero for those are inside the constraint set and they are not inside the obstable. 
function in_constraint(x, y)
    # 1. Track Boundary Check
    dx = max(abs(x) - 1.0, 0.0); dy = abs(y)
    d = sqrt(dx^2 + dy^2)
    in_track = (d <= R_outer) && (d >= R_inner)
    
    # 2. Obstacle Check
    in_obs = (x - Obs_x)^2 + (y - Obs_y)^2 <= Obs_r^2
    
    return in_track && !in_obs
end 

START_STATE  = [-2.5, -0.15, pi/4] 
TARGET_STATE_X = -1.5*ones(10)
TARGET_STATE_Y = range(-0.75, -1.75, 10) # range(start, stop; length)

if plot_test
scatter!(p,(START_STATE[1],START_STATE[2]), c=:blue, ms=4, aspect_ratio=:equal, label="Start state")
scatter!(p, TARGET_STATE_X,TARGET_STATE_Y, c=:green, ms=4, label="Target States")
display(p)
end 

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
u_limits   = [-pi/3, pi/3] # Steering Limits (rad/s)

SIGMA_SLIP = 0.1
BOUND_SLIP = 0.1
DIST_SLIP  = Truncated(Normal(0.0, SIGMA_SLIP), -BOUND_SLIP, BOUND_SLIP)

SIGMA_STEER = 0.5
BOUND_STEER = 0.5
DIST_STEER  = Truncated(Normal(0.0, SIGMA_STEER), -BOUND_STEER, BOUND_STEER)

### MDP Model 
states_3d = [[x, y, th] for y in y_grid for x in x_grid for th in th_grid]
nstates = length(states_3d)
nactions = 21
u_actions = collect(LinRange(u_limits[1], u_limits[2], nactions))

# Convert vector of vectors to matrix (3 x N) for KDTree
states_matrix = hcat(states_3d...)
tree = KDTree(states_matrix)

## Reward Distribution
function reward(s)
    return in_constraint(s[1], s[2]) ? 1.0 : 0.0
end
r_vec = [reward(s) for s in states_3d]

# ploting to check zero and one  bolean 
if plot_test
heatmap_data = [in_constraint(x, y) ? 1.0 : 0.0 for y in y_grid, x in x_grid]
c = heatmap(x_grid, y_grid, heatmap_data, 
    c=:grays, 
    aspect_ratio=:equal, 
    title="Reward Distribution for Safety Analysis",
    xlabel="X (m)", ylabel="Y (m)")
display(c)
end

# Transition Matrix 
N_SLIP_SAMPLES  = 10
N_STEER_SAMPLES = 10
TOTAL_SAMPLES   = N_SLIP_SAMPLES * N_STEER_SAMPLES
slip_samples  = rand(DIST_SLIP, N_SLIP_SAMPLES)
steer_samples = rand(DIST_STEER, N_STEER_SAMPLES)
T = [spzeros(nactions, nstates) for _ in 1:nstates]


println("Working on the transition matrix with Monte Carlo Sampling")

@time begin
    for is in 1:nstates
        s = states_3d[is]
        
        # Optimization: If state is already Unsafe (Wall/Grass), it is Absorbing.
        if !in_constraint(s[1], s[2])
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
    end
end

println("The Transition matrix with Monte Carlo Sampling is constructed.")

# checking and smoothing transition matrix 
if transition_test
    println("Negative Values checking") 
    @assert all(t->all(t.>=0.0),T) "Negative Value!"

    println("Sum is one checking") 
    for action in 1:nactions
        @assert all(i -> sum(T[i][action,:])≈1.0, axes(T,1)) "Not All row sums are approximately 1!"
    end
end 


println("Removing small values")
threshold_for_transit = 1e-4
for s in 1:nstates
	T[s]= dropzeros!(map(x->abs(x)<threshold_for_transit ? 0.0 : x, T[s]))
end


### Linear Program Optimization Primal and Dual 
#Primal 
model_p = Model(Mosek.Optimizer)
set_optimizer_attribute(model_p, "MSK_IPAR_INTPNT_BASIS", 0)

@variable(model_p, 0 <= g[1:nstates] <= 1)
@variable(model_p, h[1:nstates] >= 0)

alpha_dist = ones(nstates) ./ nstates
@objective(model_p, Min, sum(alpha_dist[s] * g[s] for s in 1:nstates))

c1 = Array{ConstraintRef}(undef, nstates, nactions)
c2 = Array{ConstraintRef}(undef, nstates, nactions)

@time begin
    for s in 1:nstates
        if r_vec[s] == 0.0

            con_g_zero = @constraint(model_p, g[s] == 0.0)
            con_h_zero = @constraint(model_p, h[s] == 0.0)

            for a in 1:nactions
                c1[s, a] = con_g_zero
                c2[s, a] = con_h_zero
            end
        else

            for a in 1:nactions
                P_row = T[s][a, :]
                
                c1[s, a] = @constraint(model_p, g[s] >= dot(P_row, g))
                c2[s, a] = @constraint(model_p, g[s] + h[s] >= 1.0 + dot(P_row, h))
            end
        end
    end
end


optimize!(model_p)
println("Dual Feasibility:")
dual_status(model_p)

println(dual_objective_value(model_p))



g_opt = value.(g)
h_opt = value.(h)


model_p[:c1] = c1
model_p[:c2] = c2

dual_matrix_c1 = value.(dual.(model_p[:c1]))
dual_matrix_c2 = value.(dual.(model_p[:c2]))

# Dual Optimization  
model_d = Model(Mosek.Optimizer)
set_optimizer_attribute(model_d, "MSK_IPAR_INTPNT_BASIS", 0)

@variable(model_d, z[1:nstates, 1:nactions] >= 0)
@variable(model_d, y[1:nstates, 1:nactions]>= 0)

@objective(model_d, Max, sum(z[s,a] for s in 1:nstates, a in 1:nactions if r_vec[s] != 0))
c1_d = Array{ConstraintRef}(undef, nstates)
c2_d = Array{ConstraintRef}(undef, nstates)

@time begin
    for s in 1:nstates
        if r_vec[s] == 0.0
            @constraint(model_d, z[s, :] .== 0.0)
            @constraint(model_d, y[s, :] .== 0.0)
            c1_d[s] = @constraint(model_d, 0.0 == 0.0)
            c2_d[s] = @constraint(model_d, 0.0 == 0.0)
        else
            c1_d[s] = @constraint(model_d, sum(z[s, :]) == 0.0)

            c2_d[s] = @constraint(model_d, sum(z[s, :]) + sum(y[s, :]) == alpha_dist[s])
        end
    end
end

@time begin
    for s in 1:nstates
        if r_vec[s] == 0.0; continue; end
        
        for a in 1:nactions
            destinations, probs = findnz(T[s][a, :])
            
            for (s_prime, prob) in zip(destinations, probs)
                if r_vec[s_prime] == 0.0; continue; end
                set_normalized_coefficient(c1_d[s_prime], z[s, a], 
                    normalized_coefficient(c1_d[s_prime], z[s, a]) - prob)

                set_normalized_coefficient(c2_d[s_prime], y[s, a], 
                    normalized_coefficient(c2_d[s_prime], y[s, a]) - prob)
            end
        end
    end
end

println("Optimization Start...")
optimize!(model_d)


z_dual=value.(z)
y_dual=value.(y)

policy = zeros(Float64, nstates, nactions)

tol = 1e-6
for s in 1:nstates
    if r_vec[s] == 0.0
        continue 
    end

    sum_z = sum(z_dual[s, :])
    sum_y = sum(y_dual[s, :])
    

    if sum_z > tol
        policy[s, :] = z_dual[s, :] ./ sum_z

    else 
        policy[s, :] = y_dual[s, :] ./ sum_y 
    end
end

println("Policy Extracted.")


# Calculate Q_safe
println("Calculate Q function based on optimal gain.")
Q_function_g = zeros(Float64, nstates, nactions)
for s in 1:nstates
    if !in_constraint(states_3d[s][1], states_3d[s][2])
        Q_function_g[s, :] .= 0.0
        continue
    end
    Q_function_g[s, :] = T[s] * g_opt
end

# Cost function 
function get_stage_cost(s)
    x, y = s[1], s[2]
    x_tgt_min = minimum(TARGET_STATE_X) - 0.1
    x_tgt_max = maximum(TARGET_STATE_X) + 0.1
    y_tgt_min = minimum(TARGET_STATE_Y)
    y_tgt_max = maximum(TARGET_STATE_Y)

    if (x_tgt_min <= x <= x_tgt_max) && (y_tgt_min <= y <= y_tgt_max)
        return 0.0  
    else 
        return 1.0 
    end
end

C_vec_base = [get_stage_cost(s) for s in states_3d]

if plot_test
    cost_grid_3d = reshape(cost_vec, num_points_th, num_points_x, num_points_y)
    cost_matrix_2d = cost_grid_3d[1, :, :]' # Transpose for (y, x) plotting

    heatmap(x_grid, y_grid, cost_matrix_2d, 
        title="Stage Cost (Dynamic Bounds)",
        xlabel="X (m)", 
        ylabel="Y (m)",
        c=:viridis,
        aspect_ratio=:equal
    )
end 



# Prepare Big Sparse Matrix for VI
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


function compute_hard_constraint_policy(alpha, name)
    println("\n   >>> Solving for '$name' (Threshold: $alpha)...")
    
    # 1. Expand State Cost to (State x Action) Dimensions
    # C_vec_base is length nstates. We repeat it for all actions.
    C_matrix_full = repeat(C_vec_base, 1, nactions) 
    C_vec_expanded = vec(C_matrix_full) # Length = nstates * nactions
    
    # 2. Mask Unsafe Actions
    valid_mask = Q_function_g .>= alpha
    
    # 3. Apply High Cost to Unsafe Actions
    C_vec_local = copy(C_vec_expanded)
    # Apply 1e6 cost where the mask is false (unsafe)
    C_vec_local[.!vec(valid_mask)] .= 1e6 

    # 4. Initialize V
    V_val = fill(1.0, nstates)
    
    # Identify target states directly from the cost vector
    # If cost is 0.0, it's a target. V(target) must be 0.0.
    target_idxs = findall(x -> x == 0.0, C_vec_base)
    V_val[target_idxs] .= 0.0
    
    gamma = 0.99
    
    # 5. Value Iteration Loop
    for k in 1:5000
        V_prev = V_val
        
        # Bellman Update: Q = C + gamma * E[V']
        # P_big is (nstates*nactions) x nstates
        Q_expected = P_big * V_prev
        Q_vec = C_vec_local .+ gamma .* Q_expected
        
        Q_mat = reshape(Q_vec, (nstates, nactions))
        
        # V(s) = min_a Q(s,a)
        V_val = minimum(Q_mat, dims=2)[:]
        
        # Enforce Target V=0 (Absorbing Goal State) to prevent drift
        V_val[target_idxs] .= 0.0
        
        if maximum(abs.(V_val - V_prev)) < 1e-2
            println("      Converged in $k iterations.")
            break
        end
    end
    
    # Final Q calculation
    Q_final_vec = C_vec_local .+ gamma .* (P_big * V_val)
    Q_final_mat = reshape(Q_final_vec, (nstates, nactions))
    
    # Generate Legal Actions List for simulation
    valid_actions_lookup = [findall(valid_mask[s, :]) for s in 1:nstates]
    
    return Q_final_mat, V_val, valid_actions_lookup
end

function transtion_matrix_policy(Q,no_sampling_condition)

    rows_r = Int[]
    cols_r = Int[]
    vals_r = Float64[]

    optimal_action_counts = zeros(Int, nstates)

    for s in 1:nstates
        if C_vec_base[s] == 0.0
            push!(rows_r, s); push!(cols_r, s); push!(vals_r, 1.0)
            continue 
        end
        if minimum(Q[s, :]) >= 1e5
            push!(rows_r, s); push!(cols_r, s); push!(vals_r, 1.0) # Self-loop
            continue 
        end
        all_optimal_actions = findall(==(minimum(Q[s, :])),Q[s, :])

        optimal_action_counts[s] = length(all_optimal_actions)

        if length(all_optimal_actions) > no_sampling_condition
            best_a = sample(all_optimal_actions,1)[1] 
        else 
            best_a = all_optimal_actions[1]
        end

        next_states, probs = findnz(T[s][best_a, :])
        
        for (ns, p) in zip(next_states, probs)
            push!(rows_r, s)   # From State
            push!(cols_r, ns)  # To State
            push!(vals_r, p)   # Probability
        end
    end

    p = sparse(rows_r, cols_r, vals_r, nstates, nstates)

    return p,optimal_action_counts
    
end


function state_disurbution(alpha_dist,P,Gamma)
    I_mat = sparse(I, nstates, nstates) # Identity Matrix for the linear solve
    A = I_mat - Gamma*P' 
    mu = A \ alpha_dist 
    mu = max.(mu, 0.0)
    return mu
end



Q_60, V_60, acts_60= compute_hard_constraint_policy(0.6, "0.60")

Q_95, V_safe_95, acts_safe_95    = compute_hard_constraint_policy(0.95,"0.95")

Q_100, V_100, acts_100    = compute_hard_constraint_policy(1.0,"1.0")



P_60, action_for_each_state_60 = transtion_matrix_policy(Q_60,5)

P_95,  action_for_each_state_95 = transtion_matrix_policy(Q_95,5)

P_100,  action_for_each_state_100 = transtion_matrix_policy(Q_100,5)

P_dual,action_for_dual= transtion_matrix_policy(policy,5)


maximum(sum(P_60, dims=2)) <= 1.0001

maximum(sum(P_95, dims=2)) <= 1.0001

maximum(sum(P_100, dims=2)) <= 1.0001

maximum(sum(P_dual, dims=2)) <= 1.0001


println("Total number of optimal policies for (alpha >= 0.6): 10^$(sum(log10.(filter(x -> x > 0, action_for_each_state_60))))")

println("Total number of optimal policies for (alpha >= 0.95): 10^$(sum(log10.(filter(x -> x > 0, action_for_each_state_95))))")

println("Total number of optimal policies for  (alpha = 1.0): 10^$(sum(log10.(filter(x -> x > 0, action_for_each_state_100))))")

println("Total number of optimal policies for dual LP policy: 10^$(sum(log10.(filter(x -> x > 0, action_for_dual))))")

# State Disturbution 
Gamma =0.99;
# alpha_dist = ones(Float64, nstates) ./ nstates

# Instead of ones(nstates)...
alpha_dist = zeros(nstates)
# Find the index of the Start State
idxs, _ = knn(tree, START_STATE, 1)
alpha_dist[idxs[1]] = 1.0

state_dis_60 = state_disurbution(alpha_dist,P_60,Gamma)
state_dis_95  = state_disurbution(alpha_dist,P_95,Gamma)
state_dis_100  = state_disurbution(alpha_dist,P_100,Gamma)

state_dual  = state_disurbution(alpha_dist,P_dual,Gamma)


# checking 
sum_60 = sum(state_dis_60)
sum_95 = sum(state_dis_95)
sum_100  = sum(state_dis_100)
sum_dual  = sum(state_dual)

println("Sum 60 Policy:\n", sum_60)
println("Sum 95 Policy: \n", sum_95)
println("Sum 100 Policy: \n", sum_100)
println("Sum Policy from Dual optimization:\n", sum_dual)


# Saving 
writedlm("state_60.csv", state_dis_60, ',')
writedlm("state_95.csv",  state_dis_95,  ',')
writedlm("state_100.csv",  state_dis_100,  ',')
writedlm("dual.csv",  state_dual,      ',')
