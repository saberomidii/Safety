using Plots
using LinearAlgebra
using Distributions
# Make sure DelimitedFiles is imported
using DelimitedFiles

"""
Implementation of the Minimum Discounted Reward Hamilton-Jacobi Formulation
for computing reachable sets, applied to a double integrator example with disturbances.
Based on the paper by Akametalu et al. (2024)
"""


# Problem parameters
const dt = 0.1               # Time step
const umax = 2.0             # Maximum control input
const discount_rate = 0.1    # Discount factor λ
# const gamma = exp(-discount_rate * dt)  # Discrete discount factor

const gamma = 0.95

# Disturbance parameters
const disturbance_mean = 0.0  # Zero mean 
const disturbance_sigma = 1.0  # Standard deviation of 1
const disturbance_dist = Normal(disturbance_mean, disturbance_sigma)
const num_disturbance_samples = 100  # Number of samples to approximate worst-case

const x1_min, x1_max = -1.0, 5.0  # State space bounds for x1
const x2_min, x2_max = -5.0, 5.0  # State space bounds for x2
const safe_x1_min, safe_x1_max = 0.0, 4.0  # Safe region bounds for x1
const safe_x2_min, safe_x2_max = -3.0, 3.0  # Safe region bounds for x2
const epsilon = 1e-5         # Convergence threshold

# Discretize the state space
const nx1 = 161              # Number of grid points for x1
const nx2 = 161              # Number of grid points for x2
const nu = 11                # Number of discrete control inputs

x1_grid = LinRange(x1_min, x1_max, nx1)
x2_grid = LinRange(x2_min, x2_max, nx2)
u_grid = LinRange(-umax, umax, nu)

# Sample disturbance values within reasonable bounds (e.g., ±3 standard deviations)
w_samples = quantile.(disturbance_dist, LinRange(0.05, 0.95, num_disturbance_samples))

# Create the target function l(x)
# For the double integrator, the target is the complement of the safe region [0,4] × [-3,3]
function target_function(x1, x2)
    # Signed distance from the safe region (negative inside, positive outside)
    dist_x1 = max(safe_x1_min - x1, x1 - safe_x1_max, 0.0)
    dist_x2 = max(safe_x2_min - x2, x2 - safe_x2_max, 0.0)
    
    # If inside safe region, return a positive value
    # If outside safe region (i.e., in target), return a negative value
    if dist_x1 == 0.0 && dist_x2 == 0.0
        return 1.0  # Inside safe region (not in target)
    else
        return -1.0  # Outside safe region (in target)
    end
end

# Compute the surface function and its clipped version
function compute_surface_function()
    l = zeros(nx1, nx2)
    for (i, x1) in enumerate(x1_grid)
        for (j, x2) in enumerate(x2_grid)
            l[i, j] = target_function(x1, x2)
        end
    end
    
    # Clip the function to [-L, L]
    L = 1.0
    h = clamp.(l .- L, -2*L, 0.0)  # h(x) = l(x) - L, clipped to [-2L, 0]
    
    return l, h
end

# Multilinear interpolation for 2D grid
function interpolate_value(x1, x2, value_grid)
    # Clamp input coordinates to be within the grid bounds
    x1_clamped = clamp(x1, x1_min, x1_max)
    x2_clamped = clamp(x2, x2_min, x2_max)
    
    # Find grid indices
    i_low = max(1, floor(Int, (x1_clamped - x1_min) / (x1_max - x1_min) * (nx1 - 1)) + 1)
    i_high = min(nx1, i_low + 1)
    j_low = max(1, floor(Int, (x2_clamped - x2_min) / (x2_max - x2_min) * (nx2 - 1)) + 1)
    j_high = min(nx2, j_low + 1)
    
    # Calculate interpolation weights
    if i_high > i_low
        x1_weight = (x1_clamped - x1_grid[i_low]) / (x1_grid[i_high] - x1_grid[i_low])
    else
        x1_weight = 0.0
    end
    
    if j_high > j_low
        x2_weight = (x2_clamped - x2_grid[j_low]) / (x2_grid[j_high] - x2_grid[j_low])
    else
        x2_weight = 0.0
    end
    
    # Clamp weights to [0, 1]
    x1_weight = clamp(x1_weight, 0.0, 1.0)
    x2_weight = clamp(x2_weight, 0.0, 1.0)
    
    # Bilinear interpolation
    val_ll = value_grid[i_low, j_low]
    val_lh = value_grid[i_low, j_high]
    val_hl = value_grid[i_high, j_low]
    val_hh = value_grid[i_high, j_high]
    
    val_l = val_ll * (1 - x2_weight) + val_lh * x2_weight
    val_h = val_hl * (1 - x2_weight) + val_hh * x2_weight
    
    return val_l * (1 - x1_weight) + val_h * x1_weight
end

# Value iteration using the Minimum Discounted Reward formulation with disturbances
function value_iteration_mdr_with_disturbance()
    # Initialize surface functions
    l, h = compute_surface_function()
    
    # Initialize value function
    U = copy(h)
    U_next = similar(U)
    
    iteration = 0
    diff = Inf
    
    while diff > epsilon
        iteration += 1
        
        for i in 1:nx1
            for j in 1:nx2
                x1 = x1_grid[i]
                x2 = x2_grid[j]
                
                # Find optimal control under worst-case disturbance
                opt_val = -Inf
                
                for u_idx in 1:nu
                    u = u_grid[u_idx]
                    
                    # Initialize min value for disturbance
                    min_disturbed_val = Inf
                    
                    # Try different disturbance values and find the worst one
                    for w in w_samples
                        # Compute next state using Euler integration with disturbance
                        next_x1 = x1 + dt * x2
                        next_x2 = x2 + dt * (u + w)
                        
                        # Interpolate value at next state
                        next_val = interpolate_value(next_x1, next_x2, U)
                        
                        # Apply discount
                        discounted_val = gamma * next_val
                        
                        # Find worst-case disturbance (minimum value)
                        min_disturbed_val = min(min_disturbed_val, discounted_val)
                    end
                    
                    # Maximize over control (robust control against worst-case disturbance)
                    opt_val = max(opt_val, min_disturbed_val)
                end
                
                # Take minimum with surface function
                U_next[i, j] = min(h[i, j], opt_val)
            end
        end
        
        # Check convergence
        diff = maximum(abs.(U_next - U))
        U .= U_next
        
        if iteration % 10 == 0
            println("Iteration $iteration, max diff: $diff")
        end
    end
    
    println("Converged after $iteration iterations")
    
    # Recover Z(x) = U(x) + L
    Z = U .+ 1.0
    
    return Z, l
end

# Plot results
function plot_results(Z, l)
    # The zero level set of Z underapproximates the reachable set
    # The safe set is the complement of the reachable set
    
    # Find the zero level set of Z
    safe_region = Z .>= 0.0
    
    p = plot(size=(800, 600), xlabel="x₁", ylabel="x₂", title="Safe Region for Double Integrator with Disturbance")
    
    # Plot the analytical safe region boundaries
    plot!(p, [safe_x1_min, safe_x1_min, safe_x1_max, safe_x1_max, safe_x1_min], 
             [safe_x2_min, safe_x2_max, safe_x2_max, safe_x2_min, safe_x2_min], 
             label="Analytical Safe Region", color=:black, linewidth=2)
    
    # Plot the computed safe region (Z ≥ 0)
    contour!(p, x1_grid, x2_grid, Z', levels=[0.0], color=:green, linewidth=3, 
             label="MDR Safe Region (Z = 0)")
    
    # Plot the computed target set level curve (l = 0)
    contour!(p, x1_grid, x2_grid, l', levels=[0.0], color=:red, linewidth=2, 
             linestyle=:dash, label="Target Set Boundary (l = 0)")
    
    savefig(p, "double_integrator_safe_region_with_disturbance.png")
    return p
end

# Function to save the computed safe region as CSV
function save_safe_region_as_csv(Z, x1_grid, x2_grid, filename="safe_region_with_disturbance.csv")
    # Open a file for writing
    open(filename, "w") do file
        # Write header
        write(file, "x1,x2,is_safe\n")
        
        # Iterate through all grid points
        for i in 1:length(x1_grid)
            for j in 1:length(x2_grid)
                x1 = x1_grid[i]
                x2 = x2_grid[j]
                
                # Check if this point is in the safe region (Z >= 0)
                is_safe = Z[i, j] >= 0.0 ? 1 : 0
                
                # Write to file
                write(file, string(x1, ",", x2, ",", is_safe, "\n"))
            end
        end
    end
    println("Safe region data saved to ", filename)
end

# Function to save just the analytical safe region as CSV
function save_analytical_safe_region_as_csv(Z, l, z_filename="Z_matrix_with_disturbance.csv", l_filename="l_matrix.csv")
    # Save Z matrix directly to CSV
    open(z_filename, "w") do file
        writedlm(file, Z, ',')
    end
    println("Z matrix saved to ", z_filename, " in ", size(Z)[1], "×", size(Z)[2], " grid format")
    
    # Save l matrix directly to CSV
    open(l_filename, "w") do file
        writedlm(file, l, ',')
    end
    println("l matrix saved to ", l_filename, " in ", size(l)[1], "×", size(l)[2], " grid format")
    
    # Save also 1D version with coordinates
    coords_z_filename = "Z_with_coords_disturbance.csv"
    open(coords_z_filename, "w") do file
        # Write header
        write(file, "x1,x2,Z_value\n")
        
        # Iterate through all grid points
        for i in 1:length(x1_grid)
            for j in 1:length(x2_grid)
                x1 = x1_grid[i]
                x2 = x2_grid[j]
                z_value = Z[i, j]
                
                # Write to file
                write(file, string(x1, ",", x2, ",", z_value, "\n"))
            end
        end
    end
    
    coords_l_filename = "l_with_coords.csv"
    open(coords_l_filename, "w") do file
        # Write header
        write(file, "x1,x2,l_value\n")
        
        # Iterate through all grid points
        for i in 1:length(x1_grid)
            for j in 1:length(x2_grid)
                x1 = x1_grid[i]
                x2 = x2_grid[j]
                l_value = l[i, j]
                
                # Write to file
                write(file, string(x1, ",", x2, ",", l_value, "\n"))
            end
        end
    end
    
    println("Coordinate versions saved to ", coords_z_filename, " and ", coords_l_filename)
end

function main()
    # Print the current working directory
    println("Current working directory: ", pwd())
    
    # First run value iteration to get Z and l with disturbance
    println("Computing value function for double integrator with Gaussian disturbance...")
    println("Using disturbance distribution: Normal(μ=", disturbance_mean, ", σ=", disturbance_sigma, ")")
    println("Sampling ", num_disturbance_samples, " disturbance values to approximate worst-case")
    
    Z, l = value_iteration_mdr_with_disturbance()
    
    # Option 1: Save the analytical safe region
    println("Saving analytical safe region data to CSV...")
    save_analytical_safe_region_as_csv(Z, l)
    
    # Option 2: Save the computed safe region
    println("Saving computed safe region data to CSV...")
    computed_file_path = joinpath(pwd(), "computed_safe_region_with_disturbance.csv")
    save_safe_region_as_csv(Z, x1_grid, x2_grid, computed_file_path)
    println("File should be at: ", computed_file_path)
    
    # Verify the file exists
    if isfile(computed_file_path)
        println("✓ File successfully created")
    else
        println("⚠ File was not created or was saved elsewhere")
    end
    
    println("Plotting results...")
    p = plot_results(Z, l)
    println("size ", size(l))
    display(p)
    
    println("Done.")
end

# Run the code
main()