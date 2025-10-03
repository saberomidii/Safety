# 2. Define the dynamics model for your specific problem
function my_dynamics_model(state::Vector{Float64}, action::Float64, disturbance::Float64)
    # If the state is not safe, return the current state (no movement)
    if !is_safe_specific(state)
        return state
    end

    x, v = state[1], state[2]
    dt = 0.1 # Time step

    # Your system equations
    x_next = x + v * dt
    v_next = v + (action + disturbance) * dt

    return [x_next, v_next]
end