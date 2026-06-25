# Packages 
using Random;Random.seed!(10)

using LinearAlgebra, SparseArrays 
using Distributions, NearestNeighbors
using JuMP, MosekTools
using Plots, LaTeXStrings
using Measures, Printf

import MathOptInterface as MOI

# Physical Model 
### Phyiscal Model 

# state grid 
xmin, xmax = 0.0, 10.0
ymin, ymax = 0.0, 5.0

dx=dy=0.1

x_grid = collect(xmin:dx:xmax)   # 101 points
y_grid = collect(ymin:dy:ymax)   # 51 points

# action rid 
umin,umax,du = -1.0,1.0,0.05
u_grid=collect(umin:du:umax)

n_steering_actions = length(u_grid)
anchor_action = n_steering_actions + 1
nactions = n_steering_actions + 1

dt,velocity = 0.1,1.

x0 =2; y0 =2.5 

# Dynamics
function zermelo_dynamics(x,y,u,wx,wy; speed=velocity)
    xn = x + dt*(speed*cos(u) + wx)
    yn = y + dt*(speed*sin(u) + wy)
    return xn,yn
end

### Sets 
# Constraint set
cxmin,cxmax = 1.0,9.0
cymin,cymax = 1.0,4.0
inside_river(x,y) = cxmin < x < cxmax && cymin < y < cymax

# Island
island_x,island_y,island_radius = 5.0,2.25,0.75
inside_island(x,y) =
    (x-island_x)^2 + (y-island_y)^2 <= island_radius^2

# Reach set
rxmin,rxmax = 8.0,8.95
rymin,rymax = 1.1,2.1
inside_reach_set(x,y) =
    rxmin <= x <= rxmax && rymin <= y <= rymax

# Constraint and avoid sets
inside_constraint_set(x,y) =
    inside_river(x,y) && !inside_island(x,y)

inside_avoid_set(x,y) =
    !inside_constraint_set(x,y)

# Shallow-water anchorage set
shallow_x,shallow_y = 4.0,3.25
shallow_a,shallow_b = 0.8,0.25

inside_shallow_water_set(x,y) =
    ((x-shallow_x)/shallow_a)^2 +
    ((y-shallow_y)/shallow_b)^2 <= 1.0 &&
    inside_constraint_set(x,y)

### MDP model 
μx,μy = -0.5,-0.75
σx,σy = 1.0,1.0
nsamples = 100

wx = rand(Normal(μx,σx),nsamples)
wy = rand(Normal(μy,σy),nsamples)

states_2d = [(x,y) for x in x_grid for y in y_grid]
nstates = length(states_2d)
tree = KDTree(hcat(collect.(states_2d)...))

T = [spzeros(nactions,nstates) for _ in 1:nstates]

@time for s in 1:nstates
    x,y = states_2d[s]

    for a in 1:nactions

        if inside_reach_set(x,y) || inside_avoid_set(x,y)

            # Reach and avoid states are absorbing
            T[s][a,s] = 1.0

        elseif a == anchor_action && inside_shallow_water_set(x,y)

            # The anchor holds inside shallow water
            T[s][a,s] = 1.0

        else

            # Steering action or unsuccessful anchoring
            u = a == anchor_action ? 0.0 : u_grid[a]
            speed = a == anchor_action ? 0.0 : velocity

            for i in 1:nsamples
                xn,yn = zermelo_dynamics(
                    x,y,u,wx[i],wy[i];
                    speed=speed
                )

                idx,_ = knn(tree,[xn,yn],1)
                T[s][a,first(idx)] += 1/nsamples
            end
        end
    end
end

println("Done building T.")
@assert all(
    isapprox(sum(T[s][a,:]),1.0;atol=1e-10)
    for s in 1:nstates, a in 1:nactions
)

println("T is valid.")

r = zeros(nstates)
c = zeros(nstates)
alpha = zeros(nstates)

for s in 1:nstates
    x,y = states_2d[s]
    r[s] = inside_reach_set(x,y)
    c[s] = inside_constraint_set(x,y)
end

idx,_ = knn(tree,[x0,y0],1)
initial_state = first(idx)
alpha[initial_state] = 1.0



# LP 
P = [spzeros(nstates,nstates) for _ in 1:nactions]
for s in 1:nstates, a in 1:nactions
    P[a][s,:]=T[s][a,:]
end

A=vcat(P...)
Iₛ=sparse(1.0I,nstates,nstates)
E=vcat([Iₛ for _ in 1:nactions]...)
r_sa=repeat(r,nactions)
c_sa=repeat(c,nactions)

function solve_safe(δ)
    model=Model(Mosek.Optimizer)
    set_silent(model)

    @variable(model,x[1:nstates*nactions]>=0)
    @variable(model,y[1:nstates*nactions]>=0)

    @objective(model,Max,r_sa⋅x)
    @constraint(model,(E-A)'*x .== 0)
    @constraint(model,E'*x+(E-A)'*y .== alpha)
    @constraint(model,c_sa⋅x>=δ)

    println("Run on: ",δ)

    optimize!(model)
    status=termination_status(model)

    status==MOI.OPTIMAL ||
        return (status=status,reward=NaN,safety=NaN,X=nothing,Y=nothing)

    xv=value.(x)
    return (
        status=status,
        reward=dot(r_sa,xv),
        safety=dot(c_sa,xv),
        X=reshape(xv,nstates,nactions),
        Y=reshape(value.(y),nstates,nactions)
    )
end

# Derive different deltas 
delta_step = 0.025 
δs = collect(0.65:delta_step:0.975)
solutions = solve_safe.(δs)

reward_values = [s.reward for s in solutions]
safety_values = [s.safety for s in solutions]
statuses      = [s.status for s in solutions]
X_values      = [s.X for s in solutions]
Y_values      = [s.Y for s in solutions]

for (δ,s) in zip(δs,solutions)
    println(
        "δ = $(round(δ,digits=2)), ",
        "reward = $(s.reward), ",
        "safety = $(s.safety), ",
        "status = $(s.status)"
    )
end

### Graph LP 
# remove NaN values
reward_values = reward_values[.!isnan.(reward_values)]
safety_values = safety_values[.!isnan.(safety_values)]


p = plot(reward_values, safety_values,
    xlabel = L"\mathbf{r}^\top \mathbf{z}",
    ylabel = L"\mathbf{c}^\top \mathbf{z}",
    xguidefont = font(14, "Computer Modern", :bold),
    yguidefont = font(14, "Computer Modern", :bold),
    linewidth = 3.5,
    aspect_ratio = 1,
    ylims = (0.65, 1.0),
    xlims = (0.1, 0.7),
    framestyle = :box,
    label = "LP solution",
    legend = :bottomleft,
    size = (640, 480),
    margin = 1mm,
    left_margin = 6mm,
    right_margin = 6mm,
    top_margin = 6mm,
    bottom_margin = 8mm
)

# Krass policy
krass_deltas = collect(0.70:0.05:0.95)
krass_ids = []

for δ in krass_deltas
    idx = findfirst(x -> isapprox(x, δ, atol=1e-6), δs)
    if idx !== nothing
        push!(krass_ids, idx)
    end
end

krass_rewards = reward_values[krass_ids]
krass_safety = safety_values[krass_ids]

println("\nKrass Policy Analysis (δ from 0.70 to 0.95):")
for (δ, r, s) in zip(krass_deltas, krass_rewards, krass_safety)
    println("δ = $(round(δ, digits=2)), reward = $(round(r, digits=4)), safety = $(round(s, digits=4))")
end 

scatter!(p, krass_rewards, krass_safety,
    marker = :circle,
    markersize = 8,
    markercolor = :red,
    markerstrokecolor = :red,
    markerstrokewidth = 1.5,
    label = "Krass policy"
)

# Avila policy
# ============================================================
# Avila-style stationary policy construction
# ============================================================

function avila_policy_from_lp(X, Y; tol=1e-10)
    ns, na = size(X)

    xmass = vec(sum(X, dims=2))
    ymass = vec(sum(Y, dims=2))

    π = zeros(ns, na)

    for s in 1:ns
        if xmass[s] > tol
            π[s, :] .= X[s, :] ./ xmass[s]
        elseif ymass[s] > tol
            π[s, :] .= Y[s, :] ./ ymass[s]
        else
            π[s, :] .= 1.0 / na
        end
    end

    return π
end


# ============================================================
# Transition matrix induced by a stationary randomized policy
# ============================================================

function policy_matrix_from_stationary_policy(π)
    Pπ = spzeros(nstates, nstates)

    for a in 1:nactions
        Pπ += spdiagm(0 => π[:, a]) * P[a]
    end

    return Pπ
end


# ============================================================
# Occupation measure induced by the Avila stationary policy
# ============================================================

function induced_occupation_measure(π)
    Pπ = policy_matrix_from_stationary_policy(π)

    model = Model(Mosek.Optimizer)
    set_silent(model)

    @variable(model, z[1:nstates] >= 0)
    @variable(model, y[1:nstates] >= 0)

    @constraint(model, (Iₛ - Pπ') * z .== 0)
    @constraint(model, z + (Iₛ - Pπ') * y .== alpha)

    @objective(model, Max, dot(r, z))

    optimize!(model)

    status = termination_status(model)

    if status != MOI.OPTIMAL
        return (
            status = status,
            reward = NaN,
            safety = NaN,
            z = nothing,
            y = nothing
        )
    end

    zv = value.(z)
    yv = value.(y)

    return (
        status = status,
        reward = dot(r, zv),
        safety = dot(c, zv),
        z = zv,
        y = yv
    )
end

# ============================================================
# Avila policy
# ============================================================

avila_deltas = collect(0.70:0.05:0.95)
avila_ids = []

for δ in avila_deltas
    idx = findfirst(x -> isapprox(x, δ, atol=1e-6), δs)

    if idx !== nothing && statuses[idx] == MOI.OPTIMAL
        push!(avila_ids, idx)
    end
end

avila_rewards = Float64[]
avila_safety = Float64[]

for idx in avila_ids
    X = X_values[idx]
    Y = Y_values[idx]

    if X !== nothing && Y !== nothing
        πA = avila_policy_from_lp(X, Y)
        valA = induced_occupation_measure(πA)

        push!(avila_rewards, valA.reward)
        push!(avila_safety, valA.safety)
    end
end

println("\nAvila Policy Analysis:")
for (δ, rA, cA) in zip(avila_deltas, avila_rewards, avila_safety)
    println(
        "δ = $(round(δ, digits=2)), ",
        "reward = $(round(rA, digits=4)), ",
        "safety = $(round(cA, digits=4))"
    )
end

scatter!(
    p,
    avila_rewards,
    avila_safety,
    marker = :diamond,
    markersize = 8,
    markercolor = :green,
    markerstrokecolor = :black,
    markerstrokewidth = 1.5,
    label = "Avila policy"
)

# Table comparing LP solution (Krass) and Avila policies
println("\n" * "="^100)
println("POLICY COMPARISON TABLE (LP Solution vs Avila Policy)")
println("="^100)

# Create header
println(
    rpad("δ", 8) * 
    rpad("LP Reward", 15) *
    rpad("LP Safety", 15) *
    rpad("Avila Reward", 15) *
    rpad("Avila Safety", 15) *
    rpad("ΔReward", 15) *
    rpad("ΔSafety", 15)
)
println("-"^100)

# Create table rows
for (i, δ) in enumerate(avila_deltas)
    lp_rew = krass_rewards[i]
    lp_safe = krass_safety[i]
    avila_rew = avila_rewards[i]
    avila_safe = avila_safety[i]
    
    delta_rew = lp_rew - avila_rew
    delta_safe = lp_safe - avila_safe
    
    println(
        rpad(string(round(δ, digits=2)), 8) *
        rpad(string(round(lp_rew, digits=4)), 15) *
        rpad(string(round(lp_safe, digits=4)), 15) *
        rpad(string(round(avila_rew, digits=4)), 15) *
        rpad(string(round(avila_safe, digits=4)), 15) *
        rpad(string(round(delta_rew, digits=4)), 15) *
        rpad(string(round(delta_safe, digits=4)), 15)
    )
end
println("="^100)
println("Note: LP Solution (Krass Policy) is the exact LP solution")
println("      ΔReward = LP Reward - Avila Reward (positive means LP performs better)")
println("      ΔSafety = LP Safety - Avila Safety (positive means LP is safer)")
println("="^100)


savefig(p, "lp_comparison.pdf")