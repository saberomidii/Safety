using Random, Distributions
using Plots
using CSV, DataFrames
using NearestNeighbors
using DelimitedFiles


# policy_path ="/Users/saber/Desktop/Safety/simulations/double_integrator/average_reward/results_primal/1.0/optimal_policy.csv"

policy_path ="/Users/saber/Desktop/Safety/Simulation/optimal_policy.csv"
o_p = CSV.read(policy_path,DataFrame, header=false)

Random.seed!(10)
optimal_policy = Matrix(o_p)

Time_step:: Int64  =500
element_type = Float64
ncols  = 1

p = zeros(element_type, Time_step, ncols)
v=  zeros(element_type, Time_step, ncols)
dt= 0.1

d_p=Normal(0.0, 1.0)

x_1_min ::Float64= -1.0
x_1_max ::Float64= 5.0
x_2_min ::Float64= -5.0
x_2_max ::Float64= 5.0

u_min ::Float64= -2.0
u_max ::Float64= 2.0

d_min::Float64=-1.0
d_max::Float64= 1.0 


num_points_action::Int64 = 41
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

p[1]= 3.0
v[1]=-1.5


no_sim = 100; 

for K=1:no_sim

    for time_index=1:Time_step-1
        next_time= time_index + 1


        d=rand(d_p,1)
        
        index_opt, distance_grid=(knn(tree,[p[time_index],v[time_index]],1))
        # println(index_opt)
        action__index=(optimal_policy[index_opt])
        # println(action__index)
        
        u=actions[action__index]
        
        # println(u)

        p[next_time]= p[time_index] + v[time_index]*dt
        v[next_time]= v[time_index] + (u[1]+d[1])*dt
    end


    # Define rectangle shape function
    rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])


    # Define constraint rectangle (change as needed)
    constraint = rectangle(4.0, 6.0, 0.0, -3.0)  # centered box for example

    # Create three subplots
    plot1 = plot(constraint, opacity=0.4, color=:red, label="Constraint Set",
                xlim=(-0.25, 4.25), ylim=(-3.5, 3.5), xlabel="Position", ylabel="Velocity",
                title="Phase Plot (p vs v)", legend=:topright)
    plot!(plot1, p, v, color=:blue, linewidth=2, label="Trajectory")

    scatter!(plot1, [2.0], [0.0], color=:green, marker=:circle, markersize=6, label="Initial State")


    plot2 = plot(p, color=:green, linewidth=2, xlabel="Time", ylabel="Position",
                title="Position Over Time", legend=false)

    plot3 = plot(v, color=:purple, linewidth=2, xlabel="Time", ylabel="Velocity",
                title="Velocity Over Time", legend=false)

    # Combine into one figure

    # assume p is 1×N
    p_to_save = vec(p')
    v_to_save = vec(v')                     # do the same for velocity

    mkpath("Safe Policy")
    CSV.write("Safe Policy/position.csv", DataFrame(position = p_to_save))
    CSV.write("Safe Policy/velocity.csv", DataFrame(velocity = v_to_save))

    final_plot=plot(plot1, plot2, plot3, layout=(3,1), size=(800, 900))

        # --- ADD THESE LINES ---
    println("Showing plot for run #$K...")
    display(final_plot) # This command forces the plot to appear on screen

    if K < no_sim
        println("Press Enter in the console to continue to the next run...")
        readline() # This pauses the script and waits for you to press Enter
    end
end

