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
num_points_state_1 = 161
num_points_state_2 = 161

x1 = collect(LinRange(x1_min, x1_max, num_points_state_1))
x2 = collect(LinRange(x2_min, x2_max, num_points_state_2))

state_2d = [(x, v) for x in x1 for v in x2]
nstates = length(state_2d)

# Control grid
u_min, u_max = -2.0, 2.0
num_points_actions = 11
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
# function s_d(x::Float64, y::Float64)
#     box_p_x = [0.0, 4.0]
#     box_p_y = [-3.0, 3.0]
#     dx = max(box_p_x[1] - x, 0.0, x - box_p_x[2])
#     dy = max(box_p_y[1] - y, 0.0, y - box_p_y[2])
#     inside_x = box_p_x[1] <= x <= box_p_x[2]
#     inside_y = box_p_y[1] <= y <= box_p_y[2]
#     if inside_x && inside_y
#         dx = minimum(abs.(x .- box_p_x))
#         dy = minimum(abs.(y .- box_p_y))
#         return -norm([dx, dy])
#     else
#         return norm([dx, dy])
#     end
# end
function s_d(x::Float64, y::Float64)
    # points x and y coordinate
    box_p_x = [0.0, 4.0]
    box_p_y = [-3.0, 3.0]

    dx = 0.0
    dy = 0.0

    if (minimum(box_p_x) <= x <= maximum(box_p_x)) &&
       (minimum(box_p_y) <= y <= maximum(box_p_y))
       
        dx = minimum([abs(x - box_p_x[1]), abs(x - box_p_x[2])])
        dy = minimum([abs(y - box_p_y[1]), abs(y - box_p_y[2])])
        m_dis = [dx, dy]
        return -norm(m_dis)

    elseif (x > maximum(box_p_x) || x < minimum(box_p_x)) &&
           (minimum(box_p_y) <= y <= maximum(box_p_y))
       
        if x > maximum(box_p_x)
            dx = abs(x - maximum(box_p_x))
        end
        if x < minimum(box_p_x)
            dx = abs(x - minimum(box_p_x))
        end
        dy = minimum([abs(y - box_p_y[1]), abs(y - box_p_y[2])])
        return norm([dx, dy])

    elseif (y > maximum(box_p_y) || y < minimum(box_p_y)) &&
           (minimum(box_p_x) <= x <= maximum(box_p_x))
       
        if y > maximum(box_p_y)
            dy = abs(y - maximum(box_p_y))
        end
        if y < minimum(box_p_y)
            dy = abs(y - minimum(box_p_y))
        end
        dx = minimum([abs(x - box_p_x[1]), abs(x - box_p_x[2])])
        return norm([dx, dy])

    else
        if x > maximum(box_p_x)
            dx = abs(x - maximum(box_p_x))
        end
        if x < minimum(box_p_x)
            dx = abs(x - minimum(box_p_x))
        end
        if y > maximum(box_p_y)
            dy = abs(y - maximum(box_p_y))
        end
        if y < minimum(box_p_y)
            dy = abs(y - minimum(box_p_y))
        end
        return norm([dx, dy])
    end
end

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
    global L = 1.0
    # h = clamp.(l .- L, -2L, 0.0)
	h = l .- L
    return l, h
end

l, value_shifted = compute_surface_function(x1, x2)
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

# ---------------------------
# Plot results
# ---------------------------
V_matrix = reshape(V, (length(x1), length(x2)))

heatmap(
    x1, x2, V_matrix',
    xlabel="x1 (position)", ylabel="x2 (velocity)",
    title="Value Function Heatmap",
    colorbar_title="V", c=:viridis
)

surface(
    x1, x2, V_matrix',
    xlabel="x1", ylabel="x2", zlabel="V",
    title="Value Function Surface"
)

V = V_matrix .+ 0*L  # Gamma :0.9900498337491681

# === Get absolute path for results folder ===
script_dir = @__DIR__                  # directory where this script is located
results_dir = joinpath(script_dir, "results")

# Create folder if it doesn't exist
if !isdir(results_dir)
    mkdir(results_dir)
    println("Created folder: $results_dir")
end

# === Save value function as CSV ===
csv_path = joinpath(results_dir, "value_function.csv")
writedlm(csv_path, V, ',')
println("Value function saved to: $csv_path")

state_csv_path = joinpath(results_dir, "state_2d.csv")
writedlm(state_csv_path, state_2d, ',')



