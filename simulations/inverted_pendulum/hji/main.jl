# Worst-case Minimum Discounted Reward for Double-Integrator
# Saber Omidi

using Distributions, LinearAlgebra, SparseArrays, Random
using NearestNeighbors, Plots, Statistics, DelimitedFiles

println("Packages imported.")

Random.seed!(10) # Sets the seed to 1234

# ---------------------------
# System and grid parameters
# ---------------------------
dt = 0.1
discount_rate = 0.1
γ = exp(-discount_rate * dt)

# State grid
x1_min, x1_max = -1.0, 5.0
x2_min, x2_max = -5.0, 5.0
num_points_state_1 = 201
num_points_state_2 = 201

x1 = collect(LinRange(x1_min, x1_max, num_points_state_1))
x2 = collect(LinRange(x2_min, x2_max, num_points_state_2))

state_2d = [(x, v) for x in x1 for v in x2]
nstates = length(state_2d)

# Control grid
u_min, u_max = -2.0, 2.0
num_points_actions = 21
actions = collect(LinRange(u_min, u_max, num_points_actions))

# Disturbances (worst-case)
sigma = 0.0
mean = 0.0
no_samples = 100
D_list = [rand(Normal(mean, sigma)) for _ in 1:no_samples]


println("Number of states: $nstates, Number of actions: $(length(actions))")

println("Disturbance list created.")
# Compute and print upper and lower bounds
d_min = minimum(D_list)
d_max = maximum(D_list)

println("Disturbance lower bound: $d_min")
println("Disturbance upper bound: $d_max")
# ---------------------------
# Signed distance function
# ---------------------------
function s_d(x::Float64, y::Float64)
    box_p_x = [0.0, 4.0]
    box_p_y = [-3.0, 3.0]
    dx = max(box_p_x[1] - x, 0.0, x - box_p_x[2])
    dy = max(box_p_y[1] - y, 0.0, y - box_p_y[2])
    inside_x = box_p_x[1] <= x <= box_p_x[2]
    inside_y = box_p_y[1] <= y <= box_p_y[2]
    if inside_x && inside_y
        dx = minimum(abs.(x .- box_p_x))
        dy = minimum(abs.(y .- box_p_y))
        return -norm([dx, dy])
    else
        return norm([dx, dy])
    end
end

# function s_d(x::Float64, y::Float64)
#     # points x and y coordinate
#     box_p_x = [0.0, 4.0]
#     box_p_y = [-3.0, 3.0]

#     dx = 0.0
#     dy = 0.0

#     if (minimum(box_p_x) <= x <= maximum(box_p_x)) &&
#        (minimum(box_p_y) <= y <= maximum(box_p_y))
       
#         dx = minimum([abs(x - box_p_x[1]), abs(x - box_p_x[2])])
#         dy = minimum([abs(y - box_p_y[1]), abs(y - box_p_y[2])])
#         m_dis = [dx, dy]
#         return -norm(m_dis)

#     elseif (x > maximum(box_p_x) || x < minimum(box_p_x)) &&
#            (minimum(box_p_y) <= y <= maximum(box_p_y))
       
#         if x > maximum(box_p_x)
#             dx = abs(x - maximum(box_p_x))
#         end
#         if x < minimum(box_p_x)
#             dx = abs(x - minimum(box_p_x))
#         end
#         dy = minimum([abs(y - box_p_y[1]), abs(y - box_p_y[2])])
#         return norm([dx, dy])

#     elseif (y > maximum(box_p_y) || y < minimum(box_p_y)) &&
#            (minimum(box_p_x) <= x <= maximum(box_p_x))
       
#         if y > maximum(box_p_y)
#             dy = abs(y - maximum(box_p_y))
#         end
#         if y < minimum(box_p_y)
#             dy = abs(y - minimum(box_p_y))
#         end
#         dx = minimum([abs(x - box_p_x[1]), abs(x - box_p_x[2])])
#         return norm([dx, dy])

#     else
#         if x > maximum(box_p_x)
#             dx = abs(x - maximum(box_p_x))
#         end
#         if x < minimum(box_p_x)
#             dx = abs(x - minimum(box_p_x))
#         end
#         if y > maximum(box_p_y)
#             dy = abs(y - maximum(box_p_y))
#         end
#         if y < minimum(box_p_y)
#             dy = abs(y - minimum(box_p_y))
#         end
#         return norm([dx, dy])
#     end
# end

# ---------------------------
# Double integrator dynamics
# ---------------------------
function di_dynamic(x::Float64, v::Float64, u::Float64, d::Float64)
    x_next = x + v*dt
    v_next = v + (u + d)*dt
    return [x_next, v_next]
end

# ---------------------------
# Surface function
# ---------------------------
function compute_surface_function(x1::Vector{Float64}, x2::Vector{Float64})
    nx1, nx2 = length(x1), length(x2)
    l = zeros(nx1, nx2)
    for (i, xi) in enumerate(x1)
        for (j, vj) in enumerate(x2)
            l[i,j] = -s_d(xi, vj)
        end
    end

    global L = maximum(abs.(l))

    # h = clamp.(l , -L, L) .- L
	h = l

    return l, h
end

l, value_shifted = compute_surface_function(x1, x2)


# surface(x1, x2, value_shifted, linewidth=2, linestyle=:dot, label="MDR under-approx")
# contour(x1, x2, value_shifted, linewidth=2, linestyle=:dot, label="MDR under-approx")

# surface(x1, x2, l, linewidth=2, linestyle=:dot, label="MDR under-approx")
# contour(x1, x2, l, linewidth=2, linestyle=:dot, label="MDR under-approx")


# contour!(x1, x2, Z', levels=[over_thr], linewidth=3, label="MDR over-approx")

V = vec(copy(value_shifted))         # flatten to vector
V_next = similar(V)

# ---------------------------
# KD-tree for nearest neighbor
# ---------------------------
states_matrix = hcat([collect(s) for s in state_2d]...)  # 2 x nstates
tree = KDTree(states_matrix)

# ---------------------------
# Value iteration: MDR worst-case
# ---------------------------
epsilon = 1e-3
diff = Inf
iteration = 0

while diff > epsilon
	global iteration
    iteration += 1

    for state_index in 1:nstates
        x, v = state_2d[state_index]
        best_over_u = -Inf

        for u in actions
            worst_over_d = Inf

            for d in D_list
                next_state = di_dynamic(x, v, u, d)
                idxs, _ = knn(tree, next_state, 1)
                j = idxs[1]

                # worst_over_d = min(worst_over_d, γ * V[j])
				worst_over_d = min(worst_over_d, γ * V[j])

            end

            best_over_u = max(best_over_u, worst_over_d)
        end

        # Surface function index using KD-tree
        idxs_surf, _ = knn(tree, [x,v], 1)
        closest_index = idxs_surf[1]

        # Minimum of surface function and worst-case value
        V_next[state_index] = min(value_shifted[closest_index], best_over_u)
    end

	global  diff 
    diff = maximum(abs.(V_next .- V))
    V .= V_next

    if iteration % 10 == 0
        println("Iteration $iteration, max diff = $diff")
    end
end

println("Value iteration converged in $iteration iterations.")



# V = V_matrix .+ 0*L

# === Get absolute path for results folder ===
script_dir = @__DIR__                  # directory where this script is located
results_dir = joinpath(script_dir, "results")

# Create folder if it doesn't exist
if !isdir(results_dir)
    mkdir(results_dir)
    println("Created folder: $results_dir")
end

# === Save value function as CSV ===
csv_path = joinpath(results_dir, "value_function_gamma_201_girds.csv")
writedlm(csv_path, V, ',')
println("Value function saved to: $csv_path")

state_csv_path = joinpath(results_dir, "state_2d.csv")
writedlm(state_csv_path, state_2d, ',')

##########################################################
# Minimum Discounted Reward (MDR) for Double-Integrator
# KNN-based Value Iteration
##########################################################

# using LinearAlgebra, Random, Distributions
# using NearestNeighbors, Plots, DelimitedFiles

# Random.seed!(123)  # reproducible

# # ---------------------------
# # 1) System parameters and grid
# # ---------------------------
# dt = 0.1
# discount_rate = 0.1
# γ = exp(-discount_rate * dt)

# # State grid
# num_points_state_1 = 501
# num_points_state_2 = 501
# x1 = collect(LinRange(-1.0, 5.0, num_points_state_1))
# x2 = collect(LinRange(-5.0, 5.0, num_points_state_2))

# # Column-major state enumeration
# state_2d = [(x, v) for v in x2, x in x1]
# nstates = length(state_2d)

# # Actions
# num_points_actions = 11
# actions = collect(LinRange(-2.0, 2.0, num_points_actions))

# println("Number of states: $nstates, Number of actions: $(length(actions))")

# # KD-tree for nearest-neighbor lookups
# states_matrix = hcat([collect(s) for s in state_2d]...)  # 2 x nstates
# tree = KDTree(states_matrix)

# # ---------------------------
# # 2) Signed distance for safe set K=[0,4]x[-3,3]
# # ---------------------------
# function s_d(x::Float64, v::Float64)
#     box_x = [0.0, 4.0]
#     box_v = [-3.0, 3.0]

#     dx = if x < box_x[1] x - box_x[1] elseif x > box_x[2] x - box_x[2] else 0.0 end
#     dv = if v < box_v[1] v - box_v[1] elseif v > box_v[2] v - box_v[2] else 0.0 end
#     outside_norm = norm([dx,dv])

#     if (box_x[1] <= x <= box_x[2]) && (box_v[1] <= v <= box_v[2])
#         dist_x = min(abs(x-box_x[1]), abs(x-box_x[2]))
#         dist_v = min(abs(v-box_v[1]), abs(v-box_v[2]))
#         return -norm([dist_x, dist_v])
#     else
#         return outside_norm
#     end
# end

# # ---------------------------
# # 3) Surface function
# # ---------------------------
# function compute_surface_function(x1, x2)
#     nx1, nx2 = length(x1), length(x2)
#     l = zeros(nx1, nx2)
#     for i in 1:nx1
#         for j in 1:nx2
#             l[i,j] = -s_d(x1[i], x2[j])
#         end
#     end
#     L = 2*maximum(abs.(l))
#     h = clamp.(l, -L, L) .- L
#     return l, h, L
# end

# # Compute and unpack
# l, value_shifted, L = compute_surface_function(x1, x2)

# # Flatten for value iteration
# V = vec(copy(value_shifted))
# V_next = similar(V)

# # ---------------------------
# # 4) Double-integrator dynamics
# # ---------------------------
# function di_dynamic(x, v, u)
#     x_next = x + v*dt
#     v_next = v + u*dt
#     return [x_next, v_next]
# end

# # ---------------------------
# # 5) MDR Value Iteration
# # ---------------------------
# epsilon = 1e-3
# diff = Inf
# iteration = 0
# max_iter = 50

# while diff > epsilon && iteration < max_iter
# 	global iteration
#     iteration += 1

#     for state_index in 1:nstates
#         x, v = state_2d[state_index]
#         best_over_u = -Inf

#         for u in actions
#             next_state = di_dynamic(x, v, u)
#             idxs, _ = knn(tree, next_state, 1)
#             val = γ * V[idxs[1]]
#             best_over_u = max(best_over_u, val)
#         end

#         # Surface function index via KNN
#         idxs_surf, _ = knn(tree, [x, v], 1)
#         V_next[state_index] = min(value_shifted[idxs_surf[1]], best_over_u)
#     end
#     global diff
#     diff = maximum(abs.(V_next .- V))
#     V .= V_next

#     if iteration % 10 == 0
#         println("Iteration $iteration, max diff = $diff")
#     end
# end

# println("Value iteration converged in $iteration iterations, final diff = $diff")

# # ---------------------------
# # 6) Plot Fig.1 style
# # ---------------------------
# V_matrix = reshape(V, length(x1), length(x2))
# Z = V_matrix .+ L

# τbar = 2.0
# under_thr = 0.0
# over_thr  = L*(1 - exp(-discount_rate*τbar))

# contour(x1, x2, Z', levels=[under_thr], linewidth=2, linestyle=:dot, label="MDR under-approx")
# contour!(x1, x2, Z', levels=[over_thr], linewidth=3, label="MDR over-approx")

# # Analytic safe set
# umax = 2.0
# xs = collect(range(0, stop=4, length=400))
# v_top = clamp.(sqrt.(2*umax*(4 .- xs)), -3, 3)
# v_bot = clamp.(-sqrt.(2*umax*(xs .- 0)), -3, 3)
# plot!(xs, v_top, color=:black, label="Analytic safe set")
# plot!(xs, v_bot, color=:black, label="")

# # Target set rectangle
# plot!([0,0,4,4,0], [-3,3,3,-3,-3], color=:red, label="Target set T")
# xlabel!("x1 (position)"); ylabel!("x2 (velocity)")
# title!("MDR Value Function Contours (Fig. 1 style)")

# # ---------------------------
# # 7) Save results
# # ---------------------------
# script_dir = @__DIR__
# results_dir = joinpath(script_dir, "results_MDR")
# isdir(results_dir) || mkdir(results_dir)

# writedlm(joinpath(results_dir, "V_matrix.csv"), V_matrix, ',')
# writedlm(joinpath(results_dir, "state_2d.csv"), state_2d, ',')

# println("Results saved to $results_dir")
