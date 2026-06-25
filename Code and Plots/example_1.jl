using Plots
using LaTeXStrings
using Measures

gr()

default(
    fontfamily = "Computer Modern",
    linewidth = 2.5,
    framestyle = :box,
    grid = false,
    guidefontsize = 14,
    tickfontsize = 12,
    legendfontsize = 10,
    titlefontsize = 14,
    size = (650, 260),
    margin = 2mm
)

# ============================================================
# Data points: [reward, safety]
# ============================================================

a2 = [0.8, 0.8]   # reward, safety
a1 = [0.0, 1.0]   # reward, safety

# ============================================================
# Plot
# ============================================================

p = plot(
    xlabel = L"\mathbf{r}^{\top}\mathbf{z}",
    ylabel = L"\mathbf{c}^{\top}\mathbf{z}",
    xlims = (-0.01, 0.81),
    ylims = (0.78, 1.01),
    xticks = 0.0:0.2:0.8,
    yticks = [0.8, 0.9, 1.0],
    framestyle = :box,
    grid = false,
    legend = :bottomleft,
    aspect_ratio = :none,
    left_margin = 4mm,
    right_margin = 2mm,
    top_margin = 2mm,
    bottom_margin = 4mm
)

plot!(
    p,
    [a1[1], a2[1]],
    [a1[2], a2[2]],
    color = :blue,
    linestyle = :dash,
    linewidth = 3.0,
    label = L"\mathrm{Combination~of~}a_1\mathrm{~and~}a_2"
)

scatter!(
    p,
    [a2[1]],
    [a2[2]],
    marker = :circle,
    markersize = 7.5,
    markercolor = :red,
    markerstrokecolor = :red,
    markerstrokewidth = 1.5,
    label = L"a_2"
)

scatter!(
    p,
    [a1[1]],
    [a1[2]],
    marker = :circle,
    markersize = 7.5,
    markercolor = :green,
    markerstrokecolor = :green,
    markerstrokewidth = 1.5,
    label = L"a_1"
)

display(p)

# ============================================================
# Save figure
# ============================================================

output_dir = "/Users/saber/Desktop/Safety/C_AVR_MDP"
mkpath(output_dir)
# savefig(p, joinpath(output_dir, "safety_optimality_tradeoff_example2_style.pdf"))
savefig(p, joinpath(output_dir, "example_1.pdf"))