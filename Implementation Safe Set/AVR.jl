module AVR

using ..Systems, ..StateSpaces
using JuMP, MosekTools, SparseArrays

export calculate_reward_vector, solve_lp, calculate_optimal_policy

"""
    calculate_reward_vector(space::StateSpace, system::AbstractSystem)

Generates a reward vector `r` where `r[s] = 1.0` if `s` is safe, and `0.0` otherwise.
"""
function calculate_reward_vector(space::StateSpace, system::AbstractSystem)
    r = zeros(Float64, space.n_states)
    for s_idx in 1:space.n_states
        if is_safe(system, space.states[s_idx])
            r[s_idx] = 1.0
        end
    end
    return r
end

"""
    solve_lp(T_SA_S, r, space, alpha_dist)

Solves the primal formulation of the average reward linear program.
Returns the objective value, g_map, h_map, and the h_opt vector.
"""
function solve_lp(
    T_SA_S::SparseMatrixCSC,
    r::Vector{Float64},
    space::StateSpace,
    alpha_dist::Vector{Float64}
)
    model = Model(Mosek.Optimizer)
    set_optimizer_attribute(model, "MSK_IPAR_INTPNT_BASIS", 0)

    @variable(model, g[1:space.n_states] >= 0)
    @variable(model, h[1:space.n_states] >= 0)

    @objective(model, Min, sum(alpha_dist[s] * g[s] for s in 1:space.n_states))

    println("Adding LP constraints...")
    for s in 1:space.n_states
        for a in 1:space.n_actions
            row_idx = (s - 1) * space.n_actions + a
            s_primes_indices, _ = findnz(T_SA_S[row_idx, :])

            if !isempty(s_primes_indices)
                @constraint(model, g[s] >= sum(g[s_prime] * T_SA_S[row_idx, s_prime] for s_prime in s_primes_indices))
                @constraint(model, g[s] + h[s] >= r[s] + sum(h[s_prime] * T_SA_S[row_idx, s_prime] for s_prime in s_primes_indices))
            end
        end
    end

    println("Solving optimization problem...")
    optimize!(model)

    if termination_status(model) in (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL)
        println("Solution found. Objective value = ", objective_value(model))
        g_opt = value.(g)
        h_opt = value.(h)
        
        g_map = reshape(g_opt, space.states[end][1] == space.states[end-1][1] ? (length(unique(s->s[1], space.states)), length(unique(s->s[2], space.states))) : (length(unique(s->s[2], space.states)), length(unique(s->s[1], space.states))))
        h_map = reshape(h_opt, space.states[end][1] == space.states[end-1][1] ? (length(unique(s->s[1], space.states)), length(unique(s->s[2], space.states))) : (length(unique(s->s[2], space.states)), length(unique(s->s[1], space.states))))

        return objective_value(model), g_map, h_map, h_opt
    else
        println("No optimal solution found. Status = ", termination_status(model))
        return nothing, nothing, nothing, nothing
    end
end

"""
    calculate_optimal_policy(T_SA_S, h_opt, space)

Extracts the optimal policy by finding the action that maximizes the
expected transient value `h` for each state.
"""
function calculate_optimal_policy(
    T_SA_S::SparseMatrixCSC,
    h_opt::Vector{Float64},
    space::StateSpace
)
    optimal_policy = zeros(Int, space.n_states)
    println("Extracting optimal policy...")
    
    for s in 1:space.n_states
        best_action_idx = 1
        max_val = -Inf
        
        for a in 1:space.n_actions
            row_idx = (s - 1) * space.n_actions + a
            s_primes_indices, _ = findnz(T_SA_S[row_idx, :])
            
            # Calculate the expected value of h for this state-action pair
            expected_h = isempty(s_primes_indices) ? 0.0 : sum(h_opt[s_prime] * T_SA_S[row_idx, s_prime] for s_prime in s_primes_indices)
            
            if expected_h > max_val
                max_val = expected_h
                best_action_idx = a
            end
        end
        optimal_policy[s] = best_action_idx
    end
    
    return optimal_policy
end

end # end module AVR