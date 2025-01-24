# ###############################################################################
# # PACKAGES
# ###############################################################################
# using Distributions
# using StatsFuns
# using JuMP
# using GLPK
# using Printf
# using Dates
# using Random
# import Gurobi
# # For 3D plotting with Makie; fully qualify calls as GLMakie.xyz
# using GLMakie

# ENV["GRB_LICENSE_FILE"] = "/Users/saber/gurobi/gurobi.lic"


# ###############################################################################
# # 1) Global Setup
# ###############################################################################
# sigma_list  = [0.2]
# center_list = [(0.0,  0.0)]

# num_points_state = 51
# positions  = collect(LinRange(-2, 2, num_points_state))
# velocities = collect(LinRange(-2, 2, num_points_state))

# num_points_action = 11
# actions = collect(LinRange(-1, 1, num_points_action))

# # 2D states: cartesian product of positions × velocities
# states = [(x, v) for x in positions for v in velocities]

# println("Number of states  = ", length(states))
# println("Number of actions = ", length(actions))


# ###############################################################################
# # 2) searchsortednearest
# ###############################################################################
# function searchsortednearest(arr::AbstractVector{<:Real}, x::Real)
#     idx = searchsortedfirst(arr, x)
#     if idx == 1
#         return 1
#     elseif idx > length(arr)
#         return length(arr)
#     else
#         if abs(arr[idx] - x) < abs(arr[idx-1] - x)
#             return idx
#         else
#             return idx - 1
#         end
#     end
# end


# ###############################################################################
# # 3) build_transition_map
# ###############################################################################
# function build_transition_map(
#     states::Vector{Tuple{Float64, Float64}},
#     actions::Vector{Float64},
#     sigma::Float64
# )
#     """
#     x_{t+1} = x + v
#     v_{t+1} = (v + u) + Normal(0, sigma)
#     Snap x_{t+1} to the nearest grid point, accumulate probability in v_{t+1}.
#     """
#     xvals = unique(sort([s[1] for s in states]))
#     vvals = unique(sort([s[2] for s in states]))
#     delta_v = abs(vvals[2] - vvals[1])

#     T = Dict{Tuple{Tuple{Float64,Float64},Float64}, Dict{Tuple{Float64,Float64},Float64}}()
#     for (x, v) in states
#         for u in actions
#             x_next_det = x + v
#             dist       = Normal(v + u, sigma)
#             local_dict = Dict{Tuple{Float64,Float64}, Float64}()

#             for vprime in vvals
#                 p_v = pdf(dist, vprime) * delta_v
#                 if p_v > 1e-15
#                     ix = searchsortednearest(xvals, x_next_det)
#                     xclosest = xvals[ix]
#                     local_dict[(xclosest, vprime)] =
#                         get(local_dict, (xclosest, vprime), 0.0) + p_v
#                 end
#             end
#             T[((x, v), u)] = local_dict
#         end
#     end
#     return T
# end


# ###############################################################################
# # 4) Reward construction
# ###############################################################################
# function is_safe(x::Float64, v::Float64)
#     return (-1 <= x <= 1) && (-1 <= v <= 1)
# end

# function build_reward(
#     states::Vector{Tuple{Float64,Float64}},
#     actions::Vector{Float64}
# )
#     R = Dict{Tuple{Tuple{Float64,Float64},Float64}, Float64}()
#     for (x,v) in states
#         for u in actions
#             R[((x,v),u)] = is_safe(x,v) ? 1.0 : 0.0
#         end
#     end
#     return R
# end

# R = build_reward(states, actions)


# ###############################################################################
# # 5) Solve One Case -> Return (X, Y, Z, c2_map, c3_map)
# ###############################################################################
# function solve_case(
#     sigma::Float64,
#     center_of_prob::Tuple{Float64,Float64}
# )

#     println("\n*********************************************")
#     @info "Solving case: σ=$sigma, α_center=$center_of_prob"
#     println("*********************************************\n")

#     @time T = build_transition_map(states, actions, sigma)

#     model = Model(Gurobi.Optimizer)
#     @variable(model, z[(s in states), a in actions] >= 0)
#     @variable(model, y[(s in states), a in actions] >= 0)

#     # c1: total measure of z is 1
#     # @constraint(model, c1, sum(z[s,a] for s in states, a in actions) == 1)

#     # c2: flow constraints
#     @constraint(model, c2[j in states],
#         sum(z[s, a] * get(T[(s, a)], j, 0.0) for s in states, a in actions)
#         == sum(z[j, a2] for a2 in actions)
#     )

#     # c3: alpha distribution
#     alpha = 1.0 / length(states)
#     @constraint(model, c3[j in states],
#         sum(z[j,a] for a in actions) +
#         sum(y[j,a] for a in actions) -
#         sum(y[s,a] * get(T[(s,a)], j, 0.0) for s in states, a in actions)
#         == alpha
#     )

#     # Objective
#     @objective(model, Max, sum(z[s,a] * R[(s,a)] for s in states, a in actions))

#     @time optimize!(model)

#     stat = termination_status(model)
#     println("Solver status: ", stat)
#     if stat == MOI.OPTIMAL
#         println("Objective value = ", objective_value(model))
#     else
#         println("No optimal solution found. status = ", stat)
#     end

#     # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     # Retrieve primal solution: z_sol
#     # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     z_sol = Dict{Tuple{Float64,Float64}, Float64}()
#     for s in states
#         # sum over actions
#         z_sol[s] = sum(value(z[s,a]) for a in actions)
#     end

#     xvals_sorted = sort(unique([s[1] for s in states]))
#     vvals_sorted = sort(unique([s[2] for s in states]))
#     Nx = length(xvals_sorted)
#     Ny = length(vvals_sorted)

#     # Build Nx×Ny matrix with occupation measure "union" \sum_a z[s,a]
#     occupation_union = zeros(Nx, Ny)
#     for i in 1:Nx
#         for j in 1:Ny
#             st = (xvals_sorted[i], vvals_sorted[j])
#             occupation_union[i,j] = z_sol[st]
#         end
#     end

#     # We will return these in "surface" format:
#     X = [xvals_sorted[i] for j in 1:Ny, i in 1:Nx]  # shape (Ny, Nx)
#     Y = [vvals_sorted[j] for j in 1:Ny, i in 1:Nx]  # shape (Ny, Nx)
#     Z = occupation_union'                           # shape (Ny, Nx)

#     # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     # Retrieve dual values for c2, c3
#     # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     primal_value_c2 = [dual(c2[j]) for j in states]
#     primal_value_c3 = [dual(c3[j]) for j in states]

#     # Reshape them into Nx×Ny
#     c2_map = reshape(primal_value_c2, Nx, Ny)
#     c3_map = reshape(primal_value_c3, Nx, Ny)

#     # We'll also return them in the same "transpose" style:
#     c2_map = c2_map'
#     c3_map = c3_map'

#     return X, Y, Z, c2_map, c3_map
# end


# ###############################################################################
# # 6) Main: Single Window, 3 Subplots with Axis & Colorbar in Sub-Grids
# ###############################################################################
# function main_3D()
#     # Create a 1×3 figure for the single case, showing:
#     #    [1,1]: occupation measure Z
#     #    [1,2]: primal_value_c2
#     #    [1,3]: primal_value_c3
#     fig = GLMakie.Figure(resolution=(1800,600))

#     s = sigma_list[1]
#     c = center_list[1]
#     X, Y, Z, c2_map, c3_map = solve_case(s, c)

#     # -- Plot 1: Occupation measure Z
#     ax1 = fig[1,1] = GLMakie.Axis3(
#         fig,  # parent figure
#         xlabel="Position (x)",
#         ylabel="Velocity (v)",
#         zlabel="z(s)",
#         title=" z(s)"
#     )
#     surf_plot1 = GLMakie.surface!(ax1, X, Y, Z, colormap=:plasma)
#     # Safe set boundary
#     safe_x = [-1, 1, 1, -1, -1]
#     safe_v = [-1, -1, 1, 1, -1]
#     safe_z = zeros(length(safe_x))
#     GLMakie.lines!(ax1, safe_x, safe_v, safe_z, color=:red, linewidth=3)

#     # -- Plot 2: primal_value_c2 (dual of flow constraint)
#     ax2 = fig[1,2] = GLMakie.Axis3(
#         fig,
#         xlabel="Position (x)",
#         ylabel="Velocity (v)",
#         zlabel="dual(c2)",
#         title="Dual regarding to c2"
#     )
#     surf_plot2 = GLMakie.surface!(ax2, X, Y, c2_map, colormap=:viridis)
#     GLMakie.lines!(ax2, safe_x, safe_v, safe_z, color=:red, linewidth=3)

#     # -- Plot 3: primal_value_c3 (dual of alpha-dist constraint)
#     ax3 = fig[1,3] = GLMakie.Axis3(
#         fig,
#         xlabel="Position (x)",
#         ylabel="Velocity (v)",
#         zlabel="dual(c3)",
#         title="Dual regarding to c3"
#     )
#     surf_plot3 = GLMakie.surface!(ax3, X, Y, c3_map, colormap=:heat)
#     GLMakie.lines!(ax3, safe_x, safe_v, safe_z, color=:red, linewidth=3)

#     # Display the figure
#     GLMakie.display(fig)
#     println("\nAll subplots shown in a single window.")
# end

# # Run main
# main_3D()

###############################################################################
# PACKAGES
###############################################################################
using Distributions
using StatsFuns
using JuMP
using GLPK
using Printf
using Dates
using Random
import Gurobi
# For 3D plotting with Makie; fully qualify calls as GLMakie.xyz
using GLMakie
using Statistics

ENV["GRB_LICENSE_FILE"] = "/Users/saber/gurobi/gurobi.lic"


###############################################################################
# 1) Global Setup
###############################################################################
sigma_list  = [0.08]
center_list = [(0.0,  0.0)]

num_points_state = 51
positions  = collect(LinRange(-2, 2, num_points_state))
velocities = collect(LinRange(-2, 2, num_points_state))


println(positions)




num_points_action = 11
actions = collect(LinRange(-1, 1, num_points_action))


# println(actions)

# # 2D states: cartesian product of positions × velocities
# states = [(x, v) for x in positions for v in velocities]

# # println(states)


# println("Number of states  = ", length(states))
# println("Number of actions = ", length(actions))


###############################################################################
# 2) searchsortednearest
###############################################################################
function searchsortednearest(arr::AbstractVector{<:Real}, x::Real)
    idx = searchsortedfirst(arr, x)
    if idx == 1
        return 1
    elseif idx > length(arr)   # out of bound case 
        return length(arr)
    else
        if abs(arr[idx] - x) < abs(arr[idx-1] - x)
            return idx
        else
            return idx - 1
        end
    end
end


###############################################################################
# 3) build_transition_map
###############################################################################
# function build_transition_map(
#     states::Vector{Tuple{Float64, Float64}},
#     actions::Vector{Float64},
#     sigma::Float64
# )
#     """
#     x_{t+1} = x + v
#     v_{t+1} = (v + u) + Normal(0, sigma)
#     Snap x_{t+1} to the nearest grid point, accumulate probability in v_{t+1}.
#     """
#     xvals = unique(sort([s[1] for s in states]))
#     vvals = unique(sort([s[2] for s in states]))
#     delta_v = abs(vvals[2] - vvals[1])

#     T = Dict{Tuple{Tuple{Float64,Float64},Float64}, Dict{Tuple{Float64,Float64},Float64}}()
#     for (x, v) in states
#         for u in actions
#             x_next_det = x + v
#             dist       = Normal(v + u, sigma)
#             local_dict = Dict{Tuple{Float64,Float64}, Float64}()

#             for vprime in vvals
#                 p_v = pdf(dist, vprime) * delta_v
#                 # println("Value for P_v",p_v)
#                 if p_v > 1e-15
#                     ix = searchsortednearest(xvals, x_next_det)
#                     xclosest = xvals[ix]
#                     local_dict[(xclosest, vprime)] =
#                         get(local_dict, (xclosest, vprime), 0.0) + p_v
#                 end
#             end
#             T[((x, v), u)] = local_dict
#         end
#     end
#     return T
# end

###############################################################################
# 4) Reward construction
###############################################################################
function is_safe(x::Float64, v::Float64)
    return (-1 <= x <= 1) && (-1 <= v <= 1)
end

# function build_reward(
#     states::Vector{Tuple{Float64,Float64}},
#     actions::Vector{Float64}
# )
#     R = Dict{Tuple{Tuple{Float64,Float64},Float64}, Float64}()
#     for (x,v) in states
#         for u in actions
#             R[((x,v),u)] = is_safe(x,v) ? 1.0 : 0.0
#         end
#     end
#     return R
# end

# R = build_reward(states, actions)


function R(s,a)

    is_safe(s[1],s[2]) ? 1.0 : 0.0

end


function  dynamics(x,v,d,u)
     if ! is_safe(x,v) 
        x_next=x
        v_next=v
     else
     x_next= x + v
     v_next = (v + u) + d
     end
     [x_next, v_next]
end


function dynamics_rand(x,v,u)
    dist = Normal(0, sigma)
    dynamics(x,v,rand(dist),u) 
end



"""
Generates all states that we will consider in the
discretized version of the problem.
"""
function gen_states()
    state_vals = Iterators.product(LinRange(-2, 2, num_points_state), LinRange(-2, 2, num_points_state))
    [sv for sv ∈ state_vals][:]
end

states :: Vector{Tuple{Float64, Float64}} = gen_states()


####

# the following replacement speeds it by about a factor of 20
using NearestNeighbors

cs2vec(s) = collect(s)

states_matrix::Matrix{Float64} = stack(cs2vec, states)
tree::KDTree = KDTree(states_matrix)

@assert first(nn(tree, cs2vec(states[30]))) == 30

####
nstates=length(states)
nactions=length(actions)


T=zeros(nstates,nactions,nstates)
nsamples=20
sigma=0.1

"""
Computes a dynamic program update for a value function (one time step)
"""

    for is ∈ 1:nstates
        s = states[is]
        ss = cs2vec(s)
        for a ∈ 1:nactions
            for i ∈ 1:nsamples
                ns = dynamics_rand(s[1],s[2], actions[a])
                nexts, dists = knn(tree, ns, 1)
                # weights = 1 ./ (dists.+0.5)
                T[is,a,first(nexts)] +=1/nsamples   
                # T[is,a,nexts] +=1/nsamples                 
            end
        end
    end

@assert (sum(T[1,1,:])) ≈ 1.0

###############################################################################
# 5) Solve One Case -> Return arrays for z(s), dual(c2) = g(s), and dual(c3).
###############################################################################
function solve_case(
    sigma::Float64,
    center_of_prob::Tuple{Float64,Float64}
)

    println("\n*********************************************")
    @info "Solving case: σ=$sigma, α_center=$center_of_prob"
    println("*********************************************\n")


    model = Model(Gurobi.Optimizer)
    @variable(model, z[(s in 1:nstates), a in 1:nactions] >= 0)
    @variable(model, y[(s in 1:nstates), a in 1:nactions] >= 0)


#     # c1: total measure of z is 1
#     # @constraint(model, c1, sum(z[s,a] for s in states, a in actions) == 1)


    # c2: flow constraints
    @constraint(model, c2[j in 1:nstates],
        sum(z[s, a] * T[s,a,j] for s in 1:nstates, a in 1:nactions)
        == sum(z[j, a2] for a2 in 1:nactions)
    )



    # c3: alpha distribution
    alpha = 1.0 / length(states)
    @constraint(model, c3[j in 1:nstates],
        sum(z[j,a] for a in 1:nactions) +
        sum(y[j,a] for a in 1:nactions) -
        sum(y[s,a] * T[s,a,j] for s in 1:nstates, a in 1:nactions)
        == alpha
    )

    # Objective
    @objective(model, Max, sum(z[s,a] * R(states[s],actions[a]) for s in 1:nstates, a in 1:nactions))

    @time optimize!(model)

    stat = termination_status(model)
    println("Solver status: ", stat)
    if stat == MOI.OPTIMAL
        println("Objective value = ", objective_value(model))
    else
        println("No optimal solution found. status = ", stat)
    end

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 1) Retrieve primal solution z(s)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    z_sol = Dict{Tuple{Float64,Float64}, Float64}()
    for s in states
        z_sol[s] = sum(value(z[s,a]) for a in actions)
    end

    xvals_sorted = sort(unique([s[1] for s in states]))
    vvals_sorted = sort(unique([s[2] for s in states]))
    Nx = length(xvals_sorted)
    Ny = length(vvals_sorted)

    # Build Nx×Ny matrix with occupation measure "union" \sum_a z[s,a]
    occupation_union = zeros(Nx, Ny)
    for i in 1:Nx
        for j in 1:Ny
            st = (xvals_sorted[i], vvals_sorted[j])
            occupation_union[i,j] = z_sol[st]
        end
    end

    # We will return these in "surface" format:
    X = [xvals_sorted[i] for j in 1:Ny, i in 1:Nx]  # shape (Ny, Nx)
    Y = [vvals_sorted[j] for j in 1:Ny, i in 1:Nx]  # shape (Ny, Nx)
    Z_occup = occupation_union'                     # shape (Ny, Nx)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 2) Retrieve dual variables for c2 (flow) and c3 (alpha-dist).
    #    Typically, c2 is "the" MDP value function g(s) in strong duality.
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    dual_c2 = [dual(c2[j]) for j in states]   # This is "g(s)" in many MDP contexts
    dual_c3 = [dual(c3[j]) for j in states]   # Another set of dual variables (often "h(s)")

    # Reshape them into Nx×Ny, then transpose
    g_map = reshape(dual_c2, Nx, Ny)'   # "the actual value function" in the usual sense
    h_map = reshape(dual_c3, Nx, Ny)'   # additional dual (like a bias or offset)

    return X, Y, Z_occup, g_map, h_map
end


###############################################################################
# 6) Main: Plot z(s), g(s), and h(s)
###############################################################################
function main_3D()
    fig = GLMakie.Figure(resolution=(1800,600))

    # Solve for the single case
    s = sigma_list[1]
    c = center_list[1]
    X, Y, Z_occup, g_map, h_map = solve_case(s, c)

    # "Safe set" boundary just for visual reference
    safe_x = [-1, 1, 1, -1, -1]
    safe_v = [-1, -1, 1, 1, -1]
    safe_z = zeros(length(safe_x))

    # --- Subplot (1,1): occupation measure z(s) ---
    ax1 = fig[1,1] = GLMakie.Axis3(
        fig,
        xlabel="Position (x)",
        ylabel="Velocity (v)",
        zlabel="z(s)",
        title="z(s)"
    )
    surf_plot1 = GLMakie.surface!(ax1, X, Y, Z_occup, colormap=:plasma)
    GLMakie.lines!(ax1, safe_x, safe_v, safe_z, color=:red, linewidth=3)

    # --- Subplot (1,2): g(s) = dual(c2) ---
    ax2 = fig[1,2] = GLMakie.Axis3(
        fig,
        xlabel="Position (x)",
        ylabel="Velocity (v)",
        zlabel="dual from constraint one",
        title="Dual 1"
    )
    surf_plot2 = GLMakie.surface!(ax2, X, Y, g_map, colormap=:viridis)
    GLMakie.lines!(ax2, safe_x, safe_v, safe_z, color=:red, linewidth=3)

    # --- Subplot (1,3): h(s) = dual(c3) ---
    #   Shown here just so you can see it. Often 'h(s)' is a bias in multi-chain problems.
    ax3 = fig[1,3] = GLMakie.Axis3(
        fig,
        xlabel="Position (x)",
        ylabel="Velocity (v)",
        zlabel="Dual from constraint two",
        title="Dual 2"
    )
    surf_plot3 = GLMakie.surface!(ax3, X, Y, h_map, colormap=:hot)
    GLMakie.lines!(ax3, safe_x, safe_v, safe_z, color=:red, linewidth=3)

    GLMakie.display(fig)
    println("\nAll subplots shown in a single 1×3 window.")
end

main_3D()




