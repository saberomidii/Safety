# MDR.jl

module MDR

using ..StateSpaces
using SparseArrays

export calculate_terminal_cost, solve_value_iteration, calculate_mdr_objective

"""
    calculate_terminal_cost(space, cost_function, L)

Calculates the terminal cost `h` by applying the provided `cost_function`
to every state in the space.
"""
function calculate_terminal_cost(space::StateSpace, cost_function::Function, L::Float64)
    h = zeros(Float64, space.n_states)
    for s_idx in 1:space.n_states
        # Apply the passed-in cost function
        dist = cost_function(space.states[s_idx])
        h[s_idx] = -dist - L
    end
    return h
end

"""
    solve_value_iteration(T_SA_S, h, space; GAMMA, MAX_ITER, TOLERANCE)

Solves for the value function using value iteration for the MDR problem.
Returns the final value function `Z` and the optimal policy.
"""
function solve_value_iteration(
    T_SA_S::SparseMatrixCSC,
    h::Vector{Float64},
    space::StateSpace;
    GAMMA::Float64,
    MAX_ITER::Int,
    TOLERANCE::Float64,
    L::Float64
)
    U = copy(h)
    policy = zeros(Int, space.n_states)
    
    println("Starting MDR value iteration...")
    for iter in 1:MAX_ITER
        U_prev = copy(U)
        
        for s in 1:space.n_states
            best_val_over_actions = -Inf
            best_action_idx = 1
            
            for a in 1:space.n_actions
                row_idx = (s - 1) * space.n_actions + a
                s_primes_indices, _ = findnz(T_SA_S[row_idx, :])
                
                min_val_over_next_states = isempty(s_primes_indices) ? -Inf : minimum(U[s_prime] for s_prime in s_primes_indices)
                
                if min_val_over_next_states > best_val_over_actions
                    best_val_over_actions = min_val_over_next_states
                    best_action_idx = a
                end
            end
            
            U[s] = min(h[s], GAMMA * best_val_over_actions)
            policy[s] = best_action_idx
        end
        
        diff = maximum(abs.(U .- U_prev))
        if iter % 100 == 0
            println("Iter: $iter, Max change: $diff")
        end
        if diff < TOLERANCE
            println("Value function converged after $iter iterations.")
            break
        end
    end
    
    Z = U .+ L
    Z_map = reshape(Z, space.states[end][1] == space.states[end-1][1] ? (length(unique(s->s[1], space.states)), length(unique(s->s[2], space.states))) : (length(unique(s->s[2], space.states)), length(unique(s->s[1], space.states))))

    return Z_map, policy
end

"""
    calculate_mdr_objective(Z_map::Matrix{Float64})

Calculates the fraction of the state space where Z >= 0.
"""
function calculate_mdr_objective(Z_map::Matrix{Float64})
    safe_state_count = count(z -> z >= 0, Z_map)
    total_states = length(Z_map)
    return safe_state_count / total_states
end

end # end module