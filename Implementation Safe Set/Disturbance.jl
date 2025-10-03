# Disturbance.jl
module DisturbanceModels

# Required packages for this module to work
using Distributions
using DelimitedFiles

# Export the names we want to be available outside of this module
export DisturbanceGenerator, generate, save_to_csv

"""
    DisturbanceGenerator

A struct to represent a disturbance generation process.
"""
struct DisturbanceGenerator
    distribution::Distribution
    min_val::Float64
    max_val::Float64
end

"""
    generate(generator::DisturbanceGenerator, nsamples::Int)

Generates a list of random disturbances based on the generator's properties.
"""
function generate(generator::DisturbanceGenerator, nsamples::Int)
    samples = rand(generator.distribution, nsamples)
    return clamp.(samples, generator.min_val, generator.max_val)
end

"""
    save_to_csv(data::AbstractVector, filepath::String)

Saves any vector of data to a specified CSV file.
"""
function save_to_csv(data::AbstractVector, filepath::String)
    writedlm(filepath, data, ',')
    println("Successfully saved data to '$(filepath)'.")
end

end # end of module