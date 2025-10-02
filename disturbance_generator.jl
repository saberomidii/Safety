using DelimitedFiles
using Distributions
using Random
Random.seed!(10)

const d_min, d_max = -1.0, 1.0

const sigma = 1.0
const mean  = 0.0
const nsamples = 100


# Create a Normal distribution object
normal_dist = Normal(0, sigma)
disturbance_list = clamp.(rand(normal_dist, nsamples),d_min,d_max)
println("No Error, Disturbance List is saved.")
writedlm("Disturbance.csv", disturbance_list, ',')
