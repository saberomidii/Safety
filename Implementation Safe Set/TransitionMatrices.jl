module TransitionMatrices

using ..Systems 
using ..StateSpaces
using SparseArrays, NearestNeighbors

export build_transition_matrix, apply_threshold

"""
    build_transition_matrix(space, system, disturbances; nsamples)

Constructs a single, large transition probability matrix of size 
(n_states * n_actions) x n_states.

This function returns the full probability matrix without any thresholding.
"""
function build_transition_matrix(
    space::StateSpace, 
    system::AbstractSystem, 
    disturbances::Vector{Float64}; 
    nsamples::Int
)
    # T_SA_S stands for Transition from (State, Action) to next State
    T_SA_S = spzeros(Float64, space.n_states * space.n_actions, space.n_states)
    n_disturbances = length(disturbances)

    println("\nBuilding transition probability matrix T...")
    Threads.@threads for s_idx in 1:space.n_states
        for a_idx in 1:space.n_actions
            # Map the (state, action) pair to a single row index
            row_idx = (s_idx - 1) * space.n_actions + a_idx
            
            # Use a dictionary to count transitions efficiently
            next_state_counts = Dict{Int, Float64}()
            current_state = space.states[s_idx]
            action = space.actions[a_idx]

            for _ in 1:nsamples
                disturbance = disturbances[rand(1:n_disturbances)]
                next_continuous_state = apply_dynamics(system, current_state, action, disturbance)
                
                idxs, _ = knn(space.knn_tree, next_continuous_state, 1)
                s_next_idx = first(idxs)
                next_state_counts[s_next_idx] = get(next_state_counts, s_next_idx, 0.0) + 1.0
            end

            # Assign probabilities to the correct row in the large sparse matrix
            for (s_next_idx, count) in next_state_counts
                T_SA_S[row_idx, s_next_idx] = count / nsamples
            end
        end
    end
    
    println("Done building T.\n")
    return T_SA_S
end

"""
    apply_threshold(T::SparseMatrixCSC, threshold::Float64)

Creates a new sparse matrix by removing entries from T that are below the threshold.
"""
function apply_threshold(T::SparseMatrixCSC, threshold::Float64)
    # Create a new matrix where small values are mapped to zero
    T_thresholded = map(x -> abs(x) < threshold ? 0.0 : x, T)
    
    # Remove the explicit zeros to ensure the matrix is perfectly sparse
    dropzeros!(T_thresholded)
    
    return T_thresholded
end

end # end module