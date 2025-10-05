# main.jl

# --- Import necessary packages ---
using DelimitedFiles, Distributions, Random, SparseArrays, LinearAlgebra, JuMP, Printf
using Plots, LaTeXStrings

# --- Include our custom modules (order matters) ---
include("Disturbance.jl")
include("Systems.jl")
include("StateSpace.jl")
include("TransitionMatrices.jl")
include("AVR.jl")
include("CostFunctions.jl")
include("MDR.jl")
include("SaveResults.jl")
include("Plotting.jl")

using .DisturbanceModels, .Systems, .StateSpaces, .TransitionMatrices, .AVR, .CostFunctions, .MDR, .SaveResults, .Plotting

# ==========================================================
# 1) CONFIGURATION
# ==========================================================
Random.seed!(10)
const RESULTS_DIR = joinpath(@__DIR__, "results")

# --- Controller Flags ---
const RUN_AVR_SOLVER = true
const RUN_MDR_SOLVER = true
const SAVE_RESULTS = true
const GENERATE_PLOT = true

# -- Disturbance Configuration --
const disturbance_type = "Normal"
const disturbance_mean = 0.0
const SIGMA = 1.0 # This is the sigma for the disturbance
const disturbance_bounds = (-1, 1)
const disturbance_nsamples = 100

# -- State and Action Space Configuration --
const x1_params = (-1.0, 5.0, 161)
const x2_params = (-5.0, 5.0, 161)
const u_params  = (-2.0, 2.0, 81)

# -- System Dynamics Configuration (for Double Integrator) --
const dt = 0.1
const safe_region = (x_lims=(0.0, 4.0), v_lims=(-3.0, 3.0))

# -- Transition Matrix Configuration --
const nsamples = 100
const threshold = 0.00

# -- MDR Solver Configuration --
const MDR_LAMBDA = 0.2;
const MDR_DT = 0.1;
const MDR_GAMMA = exp(-MDR_LAMBDA * MDR_DT)
const MDR_L_CONSTANT = sqrt((5.0 - 4.0)^2 + (5.0 - 3.0)^2)
const MDR_MAX_ITER = 100000; const MDR_TOLERANCE = 1e-9
const mdr_params = (LAMBDA = MDR_LAMBDA, DT = MDR_DT, GAMMA = MDR_GAMMA, L_CONSTANT = MDR_L_CONSTANT, MAX_ITER = MDR_MAX_ITER, TOLERANCE = MDR_TOLERANCE)

# ==========================================================
# 2) GENERATE DISTURBANCE
# ==========================================================
println("--- 1. Generating Disturbance File ---")
# --- Create a descriptive filename ---
disturbance_filename = @sprintf("Disturbance_%s_mu%.1f_sigma%.1f_b[%.1f,%.1f].csv",
    disturbance_type, disturbance_mean, SIGMA, disturbance_bounds[1], disturbance_bounds[2])

# --- Generate and save the disturbance list ---
dist_obj = Normal(disturbance_mean, SIGMA)
# Use the fully qualified name to resolve ambiguity:
dist_generator = DisturbanceModels.DisturbanceGenerator(dist_obj, disturbance_bounds[1], disturbance_bounds[2])
const disturbance_list = generate(dist_generator, disturbance_nsamples)
save_to_csv(disturbance_list, disturbance_filename)
println("Disturbance list generated and saved to '$(disturbance_filename)'.\n")

# ==========================================================
# 3) SETUP
# ==========================================================
println("--- 2. Setting up problem ---")
system = DoubleIntegrator(dt, safe_region)
space = StateSpaces.StateSpace(u_params, x1_params, x2_params)
println("System and space are ready.")

println("\n--- 3. Building Transition Matrix ---")
T_full = build_transition_matrix(space, system, disturbance_list; nsamples=nsamples)
T_sparse = apply_threshold(T_full, threshold)
println("Transition matrix built and sparsified.")

# ==========================================================
# 4) SOLVE & GATHER RESULTS
# ==========================================================
avr_data = nothing
if RUN_AVR_SOLVER
    println("\n--- Running AVR Solver ---")
    r = calculate_reward_vector(space, system)
    alpha_dist = ones(space.n_states) ./ space.n_states
    objective, g_map, h_map, h_opt = solve_lp(T_sparse, r, space, alpha_dist)
    if !isnothing(objective)
        policy = calculate_optimal_policy(T_sparse, h_opt, space)
        g1_percentage = calculate_g1_percentage(g_map)
        avr_data = (objective=objective, g_map=g_map, h_map=h_map, policy=policy, g1_percentage=g1_percentage)
        println("AVR solver finished.")
    end
end

mdr_data = nothing
if RUN_MDR_SOLVER
    println("\n--- Running MDR Solver ---")
    mdr_cost_function = make_box_cost_function(safe_region)
    h_terminal = calculate_terminal_cost(space, mdr_cost_function, MDR_L_CONSTANT)
    Z_map, policy = solve_value_iteration(T_sparse, h_terminal, space;
        GAMMA=MDR_GAMMA, MAX_ITER=MDR_MAX_ITER, TOLERANCE=MDR_TOLERANCE, L=MDR_L_CONSTANT)
    objective = calculate_mdr_objective(Z_map)
    mdr_data = (objective=objective, Z_map=Z_map, policy=policy)
    println("MDR solver finished.")
end

# ==========================================================
# 5) SAVE & PLOT RESULTS
# ==========================================================
results_dir = nothing
if SAVE_RESULTS
    disturbance_props = (type=disturbance_type, params=(mean=disturbance_mean, sigma=SIGMA),
        bounds=disturbance_bounds, nsamples=disturbance_nsamples)
    results_dir = save_analysis_results(base_dir=RESULTS_DIR, system=system, space=space,
        disturbance_props=disturbance_props, tm_nsamples=nsamples,
        avr_results=avr_data, mdr_results=mdr_data, mdr_params=mdr_params)
end

if GENERATE_PLOT && !isnothing(results_dir)
    disturbance_props_for_plot = (params=(mean=disturbance_mean, sigma=SIGMA), bounds=disturbance_bounds)
    generate_and_save_plot(filepath=joinpath(results_dir, "level_set_plot.pdf"),
        system=system, space=space, x1_params=x1_params, x2_params=x2_params,
        disturbance_props=disturbance_props_for_plot, avr_results=avr_data, mdr_results=mdr_data)
end

println("\nProgram finished.")
