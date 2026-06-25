using Random; Random.seed!(10)
using LinearAlgebra
using SparseArrays
using Distributions
using NearestNeighbors
using JuMP
using MosekTools
using Plots
using LaTeXStrings
import MathOptInterface as MOI
using Measures
using Printf

plot_sc = false  
save_fig = false 
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

if plot_sc
    θ = range(0,2π,length=300)

    p = plot(Shape([xmin,xmax,xmax,xmin],[ymin,ymin,ymax,ymax]),
        fc=:red,fa=0.15,lc=:black,label="Avoid set",
        xlims=(xmin,xmax),ylims=(ymin,ymax),aspect_ratio=:equal,
        xlabel="x",ylabel="y",legend=:topright)

    plot!(p,Shape([cxmin,cxmax,cxmax,cxmin],[cymin,cymin,cymax,cymax]),
        fc=:white,lc=:black,lw=2,label="Constraint set")

    plot!(p,Shape(shallow_x .+ shallow_a*cos.(θ),
                  shallow_y .+ shallow_b*sin.(θ)),
        fc=:lightblue,fa=0.5,lc=:blue,lw=2,label="Shallow water")

    plot!(p,Shape(island_x .+ island_radius*cos.(θ),
                  island_y .+ island_radius*sin.(θ)),
        fc=:red,fa=0.5,lc=:black,lw=2,label=false)

    plot!(p,Shape([rxmin,rxmax,rxmax,rxmin],[rymin,rymin,rymax,rymax]),
        fc=:green,fa=0.35,lc=:green,lw=2,label="Reach set")

    scatter!(p,[x0],[y0],marker=:star8,ms=9,mc=:orange,msc=:black,label="Initial state")

    annotate!(p,island_x,island_y,text("Avoid",8))
    annotate!(p,(rxmin+rxmax)/2,(rymin+rymax)/2,text("Reach",8))
    annotate!(p,shallow_x,shallow_y+shallow_b+0.15,text("Shallow",8,:blue))

    display(p)
end


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



### LP 
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
δs = collect(0.65:delta_step:0.95)
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

### Krass policy constrcution 
function krass_policy(X,Y; tol=1e-10)
    ns,na=size(X)
    xm=vec(sum(X,dims=2)); ym=vec(sum(Y,dims=2))
    front=zeros(ns,na); tail=zeros(ns,na)
    lottery=zeros(ns,2,2); rule=zeros(ns,2,na)

    for s in 1:ns
        front[s,:].=ym[s]>tol ? Y[s,:]/ym[s] : xm[s]>tol ? X[s,:]/xm[s] : 1/na
        tail[s,:] .=xm[s]>tol ? X[s,:]/xm[s] : ym[s]>tol ? Y[s,:]/ym[s] : 1/na

        lottery[s,2,:].=(0.0,1.0)

        if xm[s]>tol
            z=xm[s]+ym[s]
            lottery[s,1,:].=(ym[s]/z,xm[s]/z)
        else
            lottery[s,1,:].=(1.0,0.0)
        end

        rule[s,1,:].=front[s,:]
        rule[s,2,:].=tail[s,:]
    end

    return (
        front=front,
        tail=tail,
        lottery=lottery,
        rule=rule,
        xmass=xm,
        ymass=ym
    )
end


# Pick deltas for krass polices 
krass_deltas = [0.65,0.85,0.95]

krass_ids = [
    findfirst(d -> isapprox(d,δ; atol=1e-10),δs)
    for δ in krass_deltas
]

krass_solutions = solutions[krass_ids]

all(s.status == MOI.OPTIMAL for s in krass_solutions) ||
    error("One selected delta is infeasible.")

krass_policies = [
    krass_policy(s.X,s.Y)
    for s in krass_solutions
]


function policy_matrix(π)
    Pπ=spzeros(nstates,nstates)
    for a in 1:nactions
        Pπ+=spdiagm(0=>π[:,a])*P[a]
    end
    return Pπ
end

function krass_distribution(π; horizon=1000)
    Pf=policy_matrix(π.front)
    Pt=policy_matrix(π.tail)

    μf=copy(alpha)
    μt=zeros(nstates)

    for _ in 1:horizon
        oldf=copy(μf)
        μf=Pf'*(π.lottery[:,1,1].*oldf)
        μt=Pt'*(π.lottery[:,1,2].*oldf+μt)
    end

    return (
        front=μf,
        tail=μt,
        state=μf+μt
    )
end

krass_distributions = [
    krass_distribution(π)
    for π in krass_policies
]

state_distributions = [
    d.state for d in krass_distributions
]

state_maps = [
    reshape(d,length(y_grid),length(x_grid))
    for d in state_distributions
]

### Plot Section
gr()

default(
    fontfamily="Computer Modern",
    linewidth=2.5,
    framestyle=:box,
    grid=false,
    guidefontsize=14,
    tickfontsize=12,
    legendfontsize=10,
    titlefontsize=14,
    size=(720,360),
    margin=3mm
)

# LP plot

valid = statuses .== MOI.OPTIMAL
xplot = safety_values[valid]
yplot = reward_values[valid]

p_lp = plot(
    xplot,yplot,
    c=:blue,lw=3,
    xlabel=L"\mathbf{c}^{\top}\mathbf{z}",
    ylabel=L"\mathbf{r}^{\top}\mathbf{z}",
    label="LP solution",
    legend=:topright,
    xlims=(0.65,0.98),
    ylims=(0.0,0.75),
    xticks=([0.65,0.70,0.75,0.80,0.85,0.90,0.95], [L"0.65", L"0.70", L"0.75", L"0.80", L"0.85", L"0.90", L"0.95"])
)

plot!(
    p_lp,
    [0.65,xplot[1]],
    [yplot[1],yplot[1]],
    c=:blue,lw=3,label=false
)

for (δ,i) in zip(krass_deltas,krass_ids)
    scatter!(
        p_lp,
        [safety_values[i]],[reward_values[i]],
        mc=:black,msc=:black,ms=5,label=false
    )

    annotate!(
        p_lp,
        safety_values[i]+0.015,
        reward_values[i]+0.02,
        text(L"\delta=%$δ",9,:black)
    )
end

display(p_lp)

if save_fig
savefig(p_lp,"lp_solution_zermelo.pdf")
savefig(p_lp,"lp_solution_zermelo.png")
end 

### X and Y distribution plots
ny,nx = length(y_grid),length(x_grid)

# reshape for Plotting 
Y_maps = [
    reshape(vec(sum(s.Y,dims=2)),ny,nx)
    for s in krass_solutions
]

X_maps = [
    reshape(
        vec(sum(s.X,dims=2)),
        ny,nx
    )
    for s in krass_solutions
]

nδ = length(krass_deltas)

Ymax = maximum(maximum.(Y_maps))
Xmax = maximum(maximum.(X_maps))

nlevels = 4
Xlevels = collect(range(
    Xmax/(nlevels+1),
    Xmax*nlevels/(nlevels+1),
    length=nlevels
))

Ycmap = cgrad(
    [:white, :lavender, :plum, :mediumorchid, :purple, :black],
    [0.0, 0.05, 0.20, 0.45, 0.75, 1.0]
)

Xcolor = :red

Ymax = maximum(maximum.(Y_maps))
Yticks = collect(range(0,Ymax,length=5))
Ylabels = [@sprintf("%.2f",v) for v in Yticks]

Ycmap = cgrad([:white,:lavender,:plum,:mediumorchid,:purple,:black],[0,0.05,0.2,0.45,0.75,1])

θ = range(0,2π,length=300)

for (k,δ) in enumerate(krass_deltas)

    showleg = k == 3

    p = heatmap(
        x_grid, y_grid, Y_maps[k],
        c = Ycmap,
        clims = (0, Ymax),
        xlims = (0.9, 9.1),
        ylims = (0.9, 4.1),
        xlabel = L"x",
        ylabel = L"y",
        size = (720, 360),
        aspect_ratio = :none,
        framestyle = :box,
        grid = false,
        margin = 0mm,
        left_margin = 4mm,
        right_margin = 14mm,
        top_margin = 2mm,
        bottom_margin = 4mm,
        colorbar = true,
        legend = showleg ? :topright : false
    )

    plot!(p,
        [cxmin,cxmax,cxmax,cxmin,cxmin],
        [cymin,cymin,cymax,cymax,cymin],
        c = :black, lw = 2,
        label = showleg ? "Constraint set" : false
    )

    plot!(p,
        shallow_x .+ shallow_a*cos.(θ),
        shallow_y .+ shallow_b*sin.(θ),
        c = :cyan, ls = :dash, lw = 2,
        label = showleg ? "Shallow water" : false
    )

    plot!(p,
        Shape(
            island_x .+ island_radius*cos.(θ),
            island_y .+ island_radius*sin.(θ)
        ),
        fc = :white, lc = :black, lw = 2,
        label = false
    )

    plot!(p,
        [rxmin,rxmax,rxmax,rxmin,rxmin],
        [rymin,rymin,rymax,rymax,rymin],
        c = :green, ls = :dashdot, lw = 2,
        label = showleg ? "Reach set" : false
    )

    scatter!(p,
        [x0], [y0],
        marker = :circle, ms = 9,
        mc = :blue, msc = :black,
        label = showleg ? "Initial state" : false
    )

    annotate!(p, island_x, island_y, text("Avoid", 8, :black))
    annotate!(p, (rxmin+rxmax)/2, (rymin+rymax)/2, text("Reach", 8, :green))

    Xmass = vec(sum(krass_solutions[k].X, dims = 2))
    Xmap  = reshape(Xmass, length(y_grid), length(x_grid))

    xp   = Xmass[Xmass .> 1e-6]
    xlev = [mean(xp), maximum(xp)]

    for lev in xlev
        contour!(p,
            x_grid, y_grid, Xmap,
            levels = [lev],
            c = :darkorange,
            lw = 2.2,
            fill = false,
            label = false
        )
    end

    imax = argmax(Xmass)
    xmax_pt, ymax_pt = states_2d[imax]

    scatter!(p,
        [xmax_pt], [ymax_pt],
        marker = :circle, ms = 9,
        mc = :red, msc = :red,
        label = showleg ? L"\mathbf{z}_{max}" : false
    )

    plot!(p,
        [NaN], [NaN],
        c = :darkorange, lw = 2.2,
        label = showleg ? L"\mathbf{z}" : false
    )

    annotate!(
        p,
        ((1.18, 0.50), text(L"\mathbf{y}^{\pi}(s)", 12, :black, rotation = 90))
    )

    display(p)

    if save_fig
    δname = replace(@sprintf("%.2f",δ),"." => "_")

    savefig(p,"krass_policy_delta_$δname.pdf")
    savefig(p,"krass_policy_delta_$δname.png")
end 
end



#Print data for Table 
anchor_mask = [inside_shallow_water_set(x,y) for (x,y) in states_2d]
avoid_mask  = [inside_avoid_set(x,y)         for (x,y) in states_2d]

println(
    "δ      safety       reward       anchor       avoid        ||y||₁       zmax"
)

for (δ,sol) in zip(krass_deltas,krass_solutions)

    z = vec(sum(sol.X,dims=2))
    y = vec(sum(sol.Y,dims=2))

    safety = dot(c,z)
    reward = dot(r,z)
    anchor = sum(sol.X[anchor_mask,anchor_action])
    avoid  = sum(z[avoid_mask])
    ynorm  = sum(y)
    zmax   = maximum(z)

    @printf(
        "%.2f   %.6f   %.6f   %.6f   %.6f   %.6f   %.6f\n",
        δ,safety,reward,anchor,avoid,ynorm,zmax
    )
end