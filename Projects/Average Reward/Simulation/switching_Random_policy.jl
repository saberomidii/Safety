using Random, Distributions
using Plots
using CSV, DataFrames
using NearestNeighbors
using DelimitedFiles


cd(@__DIR__)                # guarantees pwd() == folder that holds this file

o_p = CSV.read("optimal_policy_SD_2.csv",DataFrame, header=false)

v_f = CSV.read("Value_function_SD_2.csv",DataFrame, header=false)




Random.seed!(2) 
optimal_policy = Matrix(o_p)

value_function =  Matrix(v_f)

Time_step:: Int64  =500
element_type = Float64
ncols  = 1

p = zeros(element_type, Time_step, ncols)
v=  zeros(element_type, Time_step, ncols)
dt= 0.1

d_p=Normal(0.0, 2.0)

x_1_min ::Float64= -1.0
x_1_max ::Float64= 5.0
x_2_min ::Float64= -5.0
x_2_max ::Float64= 5.0

u_min ::Float64= -2.0
u_max ::Float64= 2.0


num_points_action::Int64 = 11
num_points_state::Int64 = 161

# Number of states in each dimension
x1 = collect(LinRange(x_1_min, x_1_max, num_points_state))  # Possible x values
x2 = collect(LinRange(x_2_min, x_2_max, num_points_state))  # Possible v values

states_2d = [(x, v) for x in x1 for v in x2]



# Action discretization
actions = collect(LinRange(u_min, u_max, num_points_action))



# Build a KD-tree for snapping continuous next-states to the nearest discrete state
states_matrix = hcat([collect(s) for s in states_2d]...)  # shape: 2 x nstates
tree = KDTree(states_matrix)

p[1]= 1.0
v[1]=-1.96

random_policy_states=[]
safe_policy_states =[]

for time_index=1:Time_step - 1
    next_time= time_index  + 1

    d=rand(d_p,1)
    
    index_opt, distance_grid=(knn(tree,[p[time_index],v[time_index]],1))
    action__index=(optimal_policy[index_opt])
    
    v_value= value_function[index_opt]


    if v_value[1] >= 1 &&  (0 + 0.7 < p[time_index] <= 4 - 0.7 && -3 +0.5 < v[time_index] <= 3-0.5)

        action__index = shuffle(1:11)[1:11]
        u=actions[action__index]
        push!(random_policy_states,[p[time_index],v[time_index]])
        else   
    
        u=actions[action__index]
        push!(safe_policy_states,[p[time_index],v[time_index]])
        end
        
    p[next_time]= p[time_index] + v[time_index]*dt
    v[next_time]= v[time_index] + (u[1]+d[1])*dt
end
# # ---------------- rectangle & trajectory ----------------
# rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
# constraint = rectangle(4.0, 6.0, 0.0, -3.0)      # change if needed

# # ---------------- safe / random policy points -----------
# xs_safe = first.(safe_policy_states)
# ys_safe = last.(safe_policy_states)

# xs_rand = first.(random_policy_states)
# ys_rand = last.(random_policy_states)

# # ---------------- phase-plot (p vs v) --------------------
# plot1 = plot(constraint; opacity=0.5, color=:red, label="Constraint",
#              xlim=(-0.25, 4.25), ylim=(-3.5, 3.5),
#              xlabel="Position", ylabel="Velocity",
#              title="Phase Plot (p vs v)",
#              legend=:topright,          # place legend in top-right corner
#              legend_columns=5)          # one row, 5 items


# plot!(plot1, p, v; color=:blue, lw=2, label="Trajectory")
# scatter!(plot1, xs_safe, ys_safe; color=:green, ms=4,
#          marker=:circle,   label="Safe policy")
# scatter!(plot1, xs_rand, ys_rand; color=:orange, ms=4,
#          marker=:utriangle, label="Random policy")
# scatter!(plot1, [p[1]], [v[1]]; color=:black, ms=6,
#          marker=:star5, label="Initial state")

# # ---------------- time-series sub-plots ------------------
# plot2 = plot(p; color=:green, lw=2, xlabel="Time", ylabel="Position",
#              title="Position Over Time", legend=false)

# plot3 = plot(v; color=:purple, lw=2, xlabel="Time", ylabel="Velocity",
#              title="Velocity Over Time", legend=false)

# # ---------------- assemble the 3-row figure --------------
# plot(plot1, plot2, plot3; layout=(3,1), size=(800,900))


# # ────────── trajectory (p, v) ──────────
# p_vec = vec(p)                  # turn N×1 → N
# v_vec = vec(v)

# traj_df = DataFrame(position = p_vec,
#                     velocity = v_vec)

# CSV.write("Switching Policy/trajectory.csv", traj_df)

# # ────────── safe-policy states ─────────
# safe_df = DataFrame(
#     x = first.(safe_policy_states),  # first coordinate
#     y = last.(safe_policy_states)    # second coordinate
# )

# CSV.write("Switching Policy/safe_policy_states.csv", safe_df)

# # ────────── random-policy states ───────
# rand_df = DataFrame(
#     x = first.(random_policy_states),
#     y = last.(random_policy_states)
# )

# CSV.write("Switching Policy/random_policy_states.csv", rand_df)

