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
using GLMakie
using Statistics
using NearestNeighbors
using MeshGrid
# License Academic 
ENV["GUROBI_HOME"] = "/Library/gurobi1200/macos_universal2"
# dbceb9e4-92ab-4786-a592-f7280d754adf
# c2001602-6d04-4d8a-806a-fa87d4e8b048
###############################################################################
# 1) Grid
no_states=2
no_actions=1
d_factor=0.01

r_states_min=[-2; -2]
r_states_max=[2;2]

r_action_min=[-1]
r_action_max=[1]




states = [LinRange(r_states_min[i]:d_factor:r_states_max[i]) for i in 1:length(r_states_min)]

actions = [LinRange(r_action_min[i]:d_factor:r_action_max[i]) for i in 1:length(r_action_min)]


nsamples=20
sigma=0.1

###############################################################################
# 2) Functions
###############################################################################
function is_safe(x::Float64, v::Float64)
    return (-1 <= x <= 1) && (-1 <= v <= 1)
end

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


function Mesh_grid(states)  # two Dimension 
    state_vals = Iterators.product(states[1], states[2])
    [sv for sv ∈ state_vals][:]
end

# ###############################################################################
# # 3) Transition matrix
# ###############################################################################
grid :: Vector{Tuple{Float64, Float64}} = Mesh_grid(states)
cs2vec(s) = collect(s)
states_matrix::Matrix{Float64} = stack(cs2vec, grid)
tree::KDTree = KDTree(states_matrix)
@assert first(nn(tree, cs2vec(grid[30]))) == 30


nstates=length(grid)
nactions=length(actions)
T=zeros(nstates,nactions,nstates)
"""
Computes Transition matrix 
"""

for is ∈ 1:nstates
    s = grid[is]
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


println("Starting Optimization")


###############################################################################
# 4) Solve One Case -> Return arrays for z(s), dual(c2) = g(s), and dual(c3).
###############################################################################

println("\n*********************************************")
@info "Solving case: σ=$sigma"
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
@objective(model, Max, sum(z[s,a] * R(grid[s],actions[a]) for s in 1:nstates, a in 1:nactions))

@time optimize!(model)

stat = termination_status(model)
println("Solver status: ", stat)
if stat == MOI.OPTIMAL
    println("Objective value = ", objective_value(model))
else
    println("No optimal solution found. status = ", stat)
end

 