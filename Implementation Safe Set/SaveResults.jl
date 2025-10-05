module SaveResults

using Dates, DelimitedFiles, Printf 
using ..Systems, ..StateSpaces

export save_analysis_results

function save_analysis_results(;
    base_dir::String,
    system::AbstractSystem,
    space::StateSpace,
    disturbance_props::NamedTuple,
    tm_nsamples::Int,
    avr_results::Union{NamedTuple, Nothing} = nothing,
    mdr_results::Union{NamedTuple, Nothing} = nothing,
    mdr_params::Union{NamedTuple, Nothing} = nothing
)
    system_name = string(nameof(typeof(system)))
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    results_dir = joinpath(base_dir, system_name, timestamp)
    mkpath(results_dir)
    println("\n--- Saving results to: $(results_dir) ---")

    open(joinpath(results_dir, "summary.txt"), "w") do f
        println(f, "Safety Analysis Results")
        println(f, "========================================")
        # ... (General info is the same) ...
        println(f, "Run Timestamp: $(timestamp)")
        println(f, "System Type: $(system_name)")
        
        println(f, "\n--- System Details ---")
        for field in fieldnames(typeof(system))
            println(f, "$(field): $(getfield(system, field))")
        end

        println(f, "\n--- Discretization ---")
        println(f, "Number of States: $(space.n_states)")
        println(f, "Number of Actions: $(space.n_actions)")

        println(f, "\n--- Disturbance Properties ---")
        println(f, "Distribution Type: $(disturbance_props.type)")
        println(f, "Distribution Parameters: $(disturbance_props.params)")
        println(f, "Clamping Bounds: [$(disturbance_props.bounds[1]), $(disturbance_props.bounds[2])]")
        println(f, "Samples Generated: $(disturbance_props.nsamples)")

        println(f, "\n--- Transition Matrix ---")
        println(f, "Monte Carlo Samples per (s,a): $(tm_nsamples)")
        
        if !isnothing(avr_results)
            println(f, "\n--- Average Reward (AVR) Results ---")
            println(f, "Objective (min average reward): $(avr_results.objective)")
            @printf(f, "Percentage of states with g(x) ≈ 1: %.2f%%\n", avr_results.g1_percentage * 100)
        end

        if !isnothing(mdr_results)
            println(f, "\n--- Minimum Discounted Reward (MDR) Results ---")
            @printf(f, "Objective (Safe Area Percentage): %.2f%%\n", mdr_results.objective * 100)
            if !isnothing(mdr_params)
                println(f, "\nMDR Parameters:")
                for (param, value) in pairs(mdr_params)
                    println(f, "$(param): $(value)")
                end
            end
        end
    end

    if !isnothing(avr_results)
        writedlm(joinpath(results_dir, "G_map_AVR.csv"), avr_results.g_map, ',')
        writedlm(joinpath(results_dir, "H_map_AVR.csv"), avr_results.h_map, ',')
        writedlm(joinpath(results_dir, "policy_AVR.csv"), avr_results.policy, ',')
    end

    if !isnothing(mdr_results)
        writedlm(joinpath(results_dir, "Z_map_MDR.csv"), mdr_results.Z_map, ',')
        writedlm(joinpath(results_dir, "policy_MDR.csv"), mdr_results.policy, ',')
    end
    
    println("All results saved successfully.")
    return results_dir
end

end # end module