using JuMP, MosekTools, LinearAlgebra, MathOptInterface
using Plots, LaTeXStrings, Measures

const MOI = MathOptInterface
gr()

default(
    fontfamily="Computer Modern",
    framestyle=:box,
    grid=false,
    linewidth=2.5,
    guidefontsize=14,
    tickfontsize=12,
    legendfontsize=10,
    size=(560, 520),
    margin=1mm,
    left_margin=6mm,
    right_margin=6mm,
    top_margin=6mm,
    bottom_margin=8mm
)

# ============================================================
# MDP data
# ============================================================

ns, na = 4, 3
nZ = ns * na

P1 = Matrix{Float64}(I, ns, ns)

P2 = [
    0.00 1.00 0.00 0.00
    0.00 0.00 0.00 1.00
    0.25 0.25 0.00 0.50
    0.00 0.00 0.00 1.00
]

P3 = [
    0.00 0.00 1.00 0.00
    0.30 0.15 0.55 0.00
    0.45 0.25 0.20 0.10
    0.00 0.00 0.00 1.00
]

P = vcat(P1, P2, P3)
E = vcat([Matrix{Float64}(I, ns, ns) for _ in 1:na]...)

r_state = [0.00, 0.45, 0.82, 1.00]
c_state = [1.00, 0.85, 0.45, 0.00]

r = repeat(r_state, na)
c = repeat(c_state, na)

alpha = [0.33, 0.33, 0.34, 0.00]

# ============================================================
# LP solver
# ============================================================

function solve_lp(δ)
    model = Model(Mosek.Optimizer)
    set_silent(model)

    @variable(model, z[1:nZ] >= 0)
    @variable(model, y[1:nZ] >= 0)

    @constraint(model, (E - P)' * z .== 0)
    @constraint(model, E' * z + (E - P)' * y .== alpha)
    @constraint(model, dot(c, z) >= δ)

    @objective(model, Max, dot(r, z))
    optimize!(model)

    termination_status(model) == MOI.OPTIMAL || return nothing

    zv, yv = value.(z), value.(y)

    return (
        reward = dot(r, zv),
        safety = dot(c, zv),
        z = zv,
        y = yv
    )
end

# ============================================================
# Policy evaluation by Cesaro average
# ============================================================

function eval_policy(π; N=30000)
    Pπ = zeros(ns, ns)
    rπ = zeros(ns)
    cπ = zeros(ns)

    for s in 1:ns, a in 1:na
        idx = (a - 1) * ns + s
        Pπ[s, :] .+= π[s, a] * P[idx, :]
        rπ[s] += π[s, a] * r[idx]
        cπ[s] += π[s, a] * c[idx]
    end

    μ = copy(alpha)
    μavg = zeros(ns)

    for _ in 1:N
        μavg .+= μ
        μ = Pπ' * μ
    end

    μavg ./= N

    return (
        reward = dot(μavg, rπ),
        safety = dot(μavg, cπ)
    )
end

# ============================================================
# Deterministic stationary policies
# ============================================================

policies = [digits(k, base=na, pad=ns) .+ 1 for k in 0:(na^ns-1)]

det_reward = Float64[]
det_safety = Float64[]

for pol in policies
    π = zeros(ns, na)

    for s in 1:ns
        π[s, pol[s]] = 1.0
    end

    val = eval_policy(π)
    push!(det_reward, val.reward)
    push!(det_safety, val.safety)
end

# ============================================================
# LP frontier
# ============================================================

lp_reward = Float64[]
lp_safety = Float64[]

for δ in 0.0:0.001:1.0
    sol = solve_lp(δ)

    if sol !== nothing
        push!(lp_reward, sol.reward)
        push!(lp_safety, sol.safety)
    end
end

perm = sortperm(lp_reward)
lp_reward = lp_reward[perm]
lp_safety = lp_safety[perm]

# ============================================================
# Pick neighboring deterministic policies around δ = 0.85
# ============================================================

δline = 0.85
tol = 1e-8

left_ids = findall(det_safety .< δline - tol)
right_ids = findall(det_safety .> δline + tol)

left_safety = maximum(det_safety[left_ids])
right_safety = minimum(det_safety[right_ids])

left_same = findall(abs.(det_safety .- left_safety) .< 1e-4)
right_same = findall(abs.(det_safety .- right_safety) .< 1e-4)

iπ2 = left_same[argmax(det_reward[left_same])]
iπ1 = right_same[argmax(det_reward[right_same])]

println("π₂ = ", policies[iπ2],
        " | reward = ", det_reward[iπ2],
        " | safety = ", det_safety[iπ2])

println("π₁ = ", policies[iπ1],
        " | reward = ", det_reward[iπ1],
        " | safety = ", det_safety[iπ1])

# ============================================================
# Avila-style stationary policy from LP solution at δ = 0.85
# ============================================================

sol85 = solve_lp(δline)
zstar, ystar = sol85.z, sol85.y

πA = zeros(ns, na)

for s in 1:ns
    zsum = sum(zstar[(a - 1) * ns + s] for a in 1:na)
    ysum = sum(ystar[(a - 1) * ns + s] for a in 1:na)

    for a in 1:na
        idx = (a - 1) * ns + s

        πA[s, a] =
            zsum > tol ? zstar[idx] / zsum :
            ysum > tol ? ystar[idx] / ysum :
            1.0 / na
    end
end

valA = eval_policy(πA)

println("\nLP at δ = 0.85:")
println("reward = ", sol85.reward)
println("safety = ", sol85.safety)

println("\nAvila policy at δ = 0.85:")
println("reward = ", valA.reward)
println("safety = ", valA.safety)

# ============================================================
# Plot
# x-axis = reward, y-axis = safety
# ============================================================

p = scatter(
    det_reward,
    det_safety,
    marker=:circle,
    markersize=6.8,
    markercolor=:red,
    markerstrokecolor=:red,
    markerstrokewidth=1.2,
    alpha=0.95,
    label="DS policy",
    xlabel=L"\mathbf{r}^{\top}\mathbf{z}",
    ylabel=L"\mathbf{c}^{\top}\mathbf{z}",
    xlims=(0.0, 1.025),
    ylims=(0.0, 0.95),
    xticks=0.0:0.2:1.0,
    yticks=0.0:0.2:1.0,
    aspect_ratio=:equal,
    legend=:bottomleft,
    left_margin=8mm,
    right_margin=8mm,
    top_margin=8mm,
    bottom_margin=10mm
)

plot!(
    p,
    lp_reward,
    lp_safety,
    color=:blue,
    linewidth=3.0,
    label="LP solution"
)

hline!(
    p,
    [δline],
    color=:black,
    linestyle=:dash,
    linewidth=2.5,
    label=L"\delta = 0.85"
)

scatter!(
    p,
    [det_reward[iπ1], det_reward[iπ2]],
    [det_safety[iπ1], det_safety[iπ2]],
    marker=:circle,
    markersize=7.5,
    markercolor=:red,
    markerstrokecolor=:red,
    markerstrokewidth=1.5,
    label=false
)

scatter!(
    p,
    [valA.reward],
    [valA.safety],
    marker=:diamond,
    markersize=8,
    markercolor=:green,
    markerstrokecolor=:black,
    markerstrokewidth=1.5,
    label="Avila policy"
)

annotate!(p, det_reward[iπ2] + 0.03, det_safety[iπ2] - 0.010, text(L"\pi_1", 12, :left))
annotate!(p, det_reward[iπ1] + 0.015, det_safety[iπ1] + 0.015, text(L"\pi_2", 12, :left))
annotate!(p, valA.reward - 0.06, valA.safety - 0.02, text(L"\pi_{\mathrm{A}}", 12, :left))

display(p)

# ============================================================
# Save figure
# ============================================================

output_dir = "/Users/saber/Desktop/Safety/C_AVR_MDP"
mkpath(output_dir)

savefig(p, joinpath(output_dir, "example_2_avila_policy_switched_axes.pdf"))