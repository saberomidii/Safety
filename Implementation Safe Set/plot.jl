using DelimitedFiles
using Plots
using LaTeXStrings


plot_font = "Times New Roman"
default(fontfamily=plot_font)

cd(@__DIR__)


filePath_G = "Safety/Average Reward /Safe set for DI/results/1.0/G.csv"

filePath_Z = "Safety/Minimum discounted Reward /Safe set for DI/Z.csv"

G = readdlm(filePath_G, ',', Float64)
Z = readdlm(filePath_Z, ',', Float64)

x_bounds = range(-1, 5, length=321)
y_bounds = range(-5, 5, length=321)


box_x = [0.0, 4.0, 4.0, 0.0, 0.0] 
box_y = [-3.0, -3.0, 3.0, 3.0, -3.0]

p = contour(
    x_bounds,
    y_bounds,
    G,
    levels = [1],       
    color = :blue,
    linewidth = 2,
    colorbar = false    
)

contour!(
    p,                   
    x_bounds,
    y_bounds,
    Z,
    levels = [0],      
    color = :red,
    linewidth = 2
)

plot!([], [], color=:blue, label=L"\{ \mathbf{x} \ | \ G(\mathbf{x}) = 1 \}")
plot!([], [], color=:red, label=L"\{ \mathbf{x} \ | \ Z(\mathbf{x}) = 0 \}")
plot!(
    p,
    box_x,
    box_y,
    color = :black,
    linewidth = 2,
    label = "Constraint Set"  # Add to legend
)

annotate!(2.5, 4.35, text(L"\mathcal{D} \sim \mathcal{N}(%$(μ), %$(σ)) \ | \ \mathcal{D} \in [-1, 1]", plot_font, 14, :center))    
xlabel!(L"x_1")
ylabel!(L"x_2")


display(p)

