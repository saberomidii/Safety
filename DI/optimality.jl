# Packages 
using Random, Distributions
using NearestNeighbors
using SparseArrays
using Mosek
using JuMP
import MathOptInterface as MOI
using MosekTools  
using LinearAlgebra

using Dates
using DelimitedFiles
# Disturbance 
Random.seed!(10) # Setting the seed
μ=0.0
σ=1.0
nsamples = 100
d_list =rand(Normal(μ, σ), nsamples)
max_d = maximum(d_list)
min_d = minimum(d_list)
l_d =-1.0
up_d =1.0
d_list_bounded = clamp.(d_list,l_d,up_d)
# Transition matrix 
num_points_state_1 = 101
num_points_state_2 = 101
x1 = collect(LinRange(-1.0, 5.0, num_points_state_1))
x2 = collect(LinRange(-5.0, 5.0, num_points_state_2))
u  = collect(LinRange(-2.0, 2.0, 81))

const c_min_1 = 0.0
const c_max_1 = 4.0

const c_min_2 = -3.0
const c_max_2 =  3.0


function is_safe(x1::Float64, x2::Float64)
    # "Safe" region: x in [0,4], v in [-3,3]
    return (c_min_1 <= x1 <= c_max_1) && (c_min_2 <= x2 <= c_max_2)
end


# Continuous (noisy) dynamics for double-integrator
function di_dynamics(x1::Float64, x2::Float64, u::Float64, d::Float64)    
    if !is_safe(x1, x2)
        return (x1, x2)   # no movement if outside the safe region
    else
    dt=0.1
    x1_next = x1 + x2*dt
	x2_next = x2 + (u + d)*dt
        return (x1_next, x2_next)
    end
end


# Cartesian product of (x, v) forms the entire state space
states_2d = [(x, v) for x in x1 for v in x2]
nstates = length(states_2d)
threshold_for_transit = 0.0
# Action discretization
nactions = length(u)

# Build a KD-tree for snapping continuous next-states to the nearest discrete state
states_matrix = hcat([collect(s) for s in states_2d]...)  # shape: 3 x nstates
tree = KDTree(states_matrix)

# Initialize T array: T[s, a, s_next]
T = Vector{SparseMatrixCSC{Float64,Int64}}(undef,nstates)

for index_state in 1:nstates
	T[index_state] = spzeros(nactions,nstates)
end

# (Compute Transition Matrix)
@time begin
 for is in 1:nstates
    s = states_2d[is]   # e.g., s = (x, y, theta)
    for a in 1:nactions
        for i in 1:nsamples
                (xn, vn) = di_dynamics(s[1], s[2], u[a], d_list_bounded[i])
                # For knn, pass a 2-element Vector, not a Tuple
                idxs, dists = knn(tree, [xn, vn], 1)
                # if first(idxs) != is || j > 1
		        T[is][ a, first(idxs)] += 1.0 / nsamples
        end
    end
 end
end

println("Negative Values checking") 
@assert all(t->all(t.>=0.0),T) "Negative Value!"

println("Sum is one checking") 
for action in 1:nactions
	@assert all(i -> sum(T[i][action,:])≈1.0, axes(T,1)) "Not All row sums are approximately 1!"
end

println("Removing small values")
for s in 1:nstates
	T[s]= dropzeros!(map(x->abs(x)<threshold_for_transit ? 0.0 : x, T[s]))
end

println("Done Building T.\n")
# AVR 
# Initialize reward vector for all states
# +1 if in safe region, else 0
function reward(s::Tuple{Float64, Float64})
    return is_safe(s[1], s[2]) ? 1.0 : 0.0
end


# Initialize reward vector for all states
r = zeros(nstates)
for s in 1:nstates
    r[s] = reward(states_2d[s])
end

# LP 
# Setup model
model = Model(Mosek.Optimizer)

# TOLERANCE SETTINGS FOR A LINEAR PROGRAM
# Setting this to 0 (MSK_OFF) disables the basis identification procedure.
set_optimizer_attribute(model, "MSK_IPAR_INTPNT_BASIS", 0)
# set_time_limit_sec(model, 1800.0)

# Define variables: g[s] and h[s] for each state
@time begin
@variable(model, g[1:nstates] >= 0)
@variable(model, h[1:nstates] >= 0)
end

# Define alpha distribution (uniform for this example)
alpha_dist = ones(nstates) ./ nstates

# Objective: minimize sum_{s} alpha_s * g_s
@time begin
@objective(model, Min, sum(alpha_dist[s] * g[s] for s in 1:nstates))
end


# Precompute the indices of nonzero transitions for each (s,a)
nonzero_indices = Vector{Vector{Int}}(undef, nstates * nactions)

for s in 1:nstates
    for a in 1:nactions
        idx = (s - 1) * nactions + a
        # Find all s_prime where T[s,a,s_prime] > threshold
        nz_sprimes = findall(p -> p > threshold_for_transit, T[s][a,:])
        nonzero_indices[idx] = nz_sprimes
    end
end

@time begin
for s in 1:nstates
    for a in 1:nactions
        idx = (s - 1) * nactions + a
        s_primes = nonzero_indices[idx]
        if !isempty(s_primes)
        @constraint(model, g[s] >= sum(g[s_prime] * T[s][ a, s_prime] for s_prime in s_primes))
        end
    end
end
end

@time begin
for s in 1:nstates
    for a in 1:nactions
        idx = (s - 1) * nactions + a
        s_primes = nonzero_indices[idx]
        if !isempty(s_primes)
        @constraint(model, g[s] + h[s] >= r[s] + sum(h[s_prime] * T[s][a, s_prime] for s_prime in s_primes))
        end
    end
end
end




# print_model_summary(model)
println("\nStarting AVR LP optimization...")
# --- Start Timer for LP Solve ---
avr_lp_start_time = time()
optimize!(model)
avr_lp_solve_time = time() - avr_lp_start_time
# --------------------------------
println("AVR LP Solve Time: $(avr_lp_solve_time) seconds")

# Check solution status
stat = termination_status(model)
println("Solver status: $stat")

println("Solution found with status: $stat")
println("Objective value = ", objective_value(model))
            
# Retrieve the g and h values (optimal or best found)
g_opt = value.(g)
h_opt = value.(h)

# Check for NaN or Inf values and replace them
# g_opt = replace(g_opt, NaN => 0.0, Inf => 1.0, -Inf => 0.0)
# h_opt = replace(h_opt, NaN => 0.0, Inf => 1.0, -Inf => 0.0)
# Reshape for visualization
gain = reshape(g_opt, num_points_state_1, num_points_state_2)
bias = reshape(h_opt, num_points_state_1, num_points_state_2)
# Extract optimal policy (for each state, find action that maximizes the value)
optimal_policy = zeros(Int, nstates)
@time optimal_policy = [argmax(r[s] .+ (T[s] * h_opt)) for s in 1:nstates]
policy= reshape(optimal_policy,num_points_state_1,num_points_state_2)

println("--- AVR is done")

## State action space 
Q_g_star = zeros(nstates, nactions)
for s in 1:nstates
    Q_g_star[s, :] = r[s] .+ (T[s] * g_opt)
end

Q_h = zeros(nstates, nactions)
for s in 1:nstates
        Q_h[s, :] = r[s] .+ (T[s] * h_opt)
end
println("Q functions are calculated.")


#Constrained Value iteration
lambda=0.05
dt=0.1
gamma = exp(-lambda*dt)       
max_iter = 2000    
epsilon = 1e-6    
safety_tol = 1e-4   

# Quadratic Cost 
Q_pos = 1.0
Q_vel = 1.0
R_u   = 0.1

# Target state 
target_x = 2.0
target_v = 0.0

# --- Initialization ---
V = zeros(nstates)              
V_next = zeros(nstates)       
optimal_policy_idx = zeros(Int, nstates) 

println("Starting Constrained Value Iteration...")
for iter in 1:max_iter
    max_diff = 0.0
    
    for s in 1:nstates
        best_val = 1e8 
        best_a = 0

        (pos, vel) = states_2d[s]
        for a in 1:nactions

            if Q_g_star[s, a] < (1.0 - safety_tol)
                continue 
            end
            # 2. Calculate Immediate Cost (Quadratic)
            action_val = u[a]
            stage_cost = Q_pos * (pos - target_x)^2 + 
                            Q_vel * (vel - target_v)^2 + 
                            R_u * (action_val)^2

            expected_future_val = dot(T[s][a, :], V)
            
            # 4. Bellman Equation
            q_value = stage_cost + (gamma * expected_future_val)
            
            # 5. Minimization
            if q_value < best_val
                best_val = q_value
                best_a = a
            end
        end
        
        # Update buffer
        V_next[s] = best_val
        
        # Store optimal action index if a safe one was found
        if best_a > 0
            optimal_policy_idx[s] = best_a
        end
        
        # Track max change for convergence
        diff = abs(V_next[s] - V[s])
        if diff > max_diff
            max_diff = diff
        end
    end
    
    # Synchronous update of Value Function
    V .= V_next
    
    # Check Convergence
    if max_diff < epsilon
        println("Converged at iteration $iter with max diff: $max_diff")
        break
    end
    
    if iter % 50 == 0
        println("Iteration $iter, max_diff: $max_diff")
    end
end

# Convert optimal indices to actual control values (u)
# If index is 0 (no safe action), we assign NaN
optimal_policy_u = [idx > 0 ? u[idx] : NaN for idx in optimal_policy_idx]

println("Optimization complete.")
println("Safe Optimal Policy stored in 'optimal_policy_u'.")


writedlm(joinpath(@__DIR__, "value_function.csv"), V, ',')
writedlm(joinpath(@__DIR__, "optimal_policy.csv"), optimal_policy_u, ',')




# Policy iteration test 
println("\n--- Starting Constrained Policy Iteration ---")

max_pi_iters = 50       
max_eval_iters = 1000   
eval_epsilon = 1e-2     
safety_tol = 1e-4      

# --- Initialization ---
# 1. Initialize Value Function and Policy
V_pi = zeros(nstates)
policy_idx = zeros(Int, nstates) 

println("Initializing with a default safe policy...")
for s in 1:nstates
    for a in 1:nactions
        if Q_g_star[s, a] >= (1.0 - safety_tol)
            policy_idx[s] = a
            break
        end
    end
    if policy_idx[s] == 0
        V_pi[s] = 1e8 
    end
end

# --- Policy Iteration Loop ---
for pi_iter in 1:max_pi_iters
    
    for eval_k in 1:max_eval_iters
        max_eval_diff = 0.0
        
        for s in 1:nstates
            current_action_idx = policy_idx[s]
            
            if current_action_idx == 0
                continue
            end
            
            # A. Immediate Cost
            (pos, vel) = states_2d[s]
            u_val = u[current_action_idx]
            cost = Q_pos*(pos - target_x)^2 + Q_vel*(vel - target_v)^2 + R_u*(u_val)^2
            
            # B. Expected Future Value (using current policy's transition)
            # T[s] is (nactions x nstates), we select the row for current_action_idx
            expected_val = dot(T[s][current_action_idx, :], V_pi)
            
            # C. Update V
            v_new = cost + gamma * expected_val
            
            diff = abs(v_new - V_pi[s])
            if diff > max_eval_diff
                max_eval_diff = diff
            end
            
            V_pi[s] = v_new
        end
        
        if max_eval_diff < eval_epsilon
            # Converged for this specific policy
            break
        end
    end

    # STEP 2: POLICY IMPROVEMENT
    policy_stable = true
    changes = 0
    for s in 1:nstates
        old_action = policy_idx[s]
        
        # Skip if state is unsafe
        if old_action == 0 
            continue 
        end
        
        best_val = Inf
        best_action = old_action
        
        (pos, vel) = states_2d[s]
        
        # Search ALL actions to find if a better one exists
        for a in 1:nactions
            
            # --- SAFETY CONSTRAINT ---
            if Q_g_star[s, a] < (1.0 - safety_tol)
                continue
            end
            
            # Calculate Q(s,a)
            u_val = u[a]
            cost = Q_pos*(pos - target_x)^2 + Q_vel*(vel - target_v)^2 + R_u*(u_val)^2
            expected_val = dot(T[s][a, :], V_pi)
            
            q_val = cost + gamma * expected_val
            
            if q_val < best_val
                best_val = q_val
                best_action = a
            end
        end
        
        # Check if policy changed
        if best_action != old_action
            policy_idx[s] = best_action
            policy_stable = false
            changes += 1
        end
    end
    
    println("PI Iteration $pi_iter: Policy changes = $changes")
    
    if policy_stable
        println("Policy Iteration Converged at iteration $(pi_iter)!")
        break
    end
end


# --- Finalize Data ---
optimal_policy_u_pi = [idx > 0 ? u[idx] : NaN for idx in policy_idx]

println("Policy Iteration Complete.")

# --- Save Results ---
writedlm(joinpath(@__DIR__, "value_function_pi.csv"), V_pi, ',')
writedlm(joinpath(@__DIR__, "optimal_policy_pi.csv"), optimal_policy_u_pi, ',')
println("PI Data saved.")