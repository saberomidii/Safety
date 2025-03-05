using DelimitedFiles, GLMakie

# Define the folder name
# Correct way to define the folder path
data_folder = "/Users/saber/Desktop/safety/Projects/Average Reward/min_d_reward_general_safety/results_primal/1.0/"


println(isfile(joinpath(data_folder, "X_grid.csv")))

# Alternative: Using raw string (prevents escape issues)
# data_folder = raw"C:\Users\saber\Desktop\safety\Projects\Average Reward\Inverted Pendulum\Data_driven_paper\0.25"



# Read CSV files from the folder
X = readdlm(joinpath(data_folder, "X_grid.csv"), ',')
Y = readdlm(joinpath(data_folder, "Y_grid.csv"), ',')
G_map = readdlm(joinpath(data_folder, "G_map.csv"), ',')
H_map = readdlm(joinpath(data_folder, "H_map.csv"), ',')

##########################################################
# 5) Plot in 3D (z, g, h) side-by-side (No Colorbars)
#    Then save results to a new directory "results/"
##########################################################
function main_3D(X, Y, G_map, H_map)
    # "Safe set" boundary lines (for reference)
    safe_x = [0, 4, 4, 0, 0]
    safe_v = [-3, -3, 3, 3, -3]
    safe_z = zeros(length(safe_x))
 
    # Create a figure with 3 subplots (no colorbars)
    fig = Figure(resolution=(1800, 600))


    # (1) Dual w.r.t. flow constraint: G
    ax2 = Axis3(fig[1,1],
        title  = "G",
        xlabel = "x",
        ylabel = "v",
        zlabel = "Value"
    )
    surface!(ax2, X, Y, G_map', colormap=:viridis)
    lines!(ax2, safe_x, safe_v, safe_z, color=:red, linewidth=3)

    # (3) Dual w.r.t. alpha-dist constraint: H_map
    ax3 = Axis3(fig[1,2],
        title  = "H",
        xlabel = "x",
        ylabel = "v",
        zlabel = "Value"
    )
    surface!(ax3, X, Y, H_map', colormap=:hot)
    lines!(ax3, safe_x, safe_v, safe_z, color=:red, linewidth=3)

    display(fig)
    println("\nAll subplots shown in one row (no colorbars).")

    # -------------------------------------------------------
    # Create a "results" directory, save the figure there.
    # -------------------------------------------------------
    mkpath("results")  # Create the directory if it doesn't exist
    GLMakie.save("results/plot.png", fig)  # Save as a .png image
    println("Figure saved to results/plot.png")

    return fig
end

##########################################################
# 6) Run
##########################################################
# Run the function with loaded data
main_3D(X, Y, G_map, H_map)
