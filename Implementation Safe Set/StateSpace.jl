# StateSpace.jl
module StateSpaces

using NearestNeighbors

export StateSpace

"""
    StateSpace

Holds the discretized state space, action space, and a k-NN tree for lookups.
"""
struct StateSpace
    states::Vector{Vector{Float64}}
    actions::Vector{Float64}
    n_states::Int
    n_actions::Int
    knn_tree::KDTree

    # This is a flexible constructor that can build a state space of any dimension.
    function StateSpace(u_range, state_ranges...)
        # Create action space from the first argument
        actions = collect(LinRange(u_range...))
        n_actions = length(actions)

        # Create state space dimensions from the rest of the arguments
        state_dims = [collect(LinRange(r...)) for r in state_ranges]
        
        # Create the Cartesian product of all state dimensions
        states_as_tuples = Iterators.product(state_dims...)
        
        # Convert states from tuples to vectors. The `vec()` function flattens
        # the resulting matrix of states into the required vector of states.
        states = vec([collect(s) for s in states_as_tuples])
        n_states = length(states)
        
        # Build the KD-tree. The hcat function needs a splatted array of vectors.
        states_matrix = hcat(states...) 
        knn_tree = KDTree(states_matrix)
        
        new(states, actions, n_states, n_actions, knn_tree)
    end
end

end # end module