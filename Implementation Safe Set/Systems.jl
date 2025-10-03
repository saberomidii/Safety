module Systems

export AbstractSystem, is_safe, apply_dynamics
export DoubleIntegrator, InvertedPendulum, DubinsCar

# Define an abstract type for any dynamic system
abstract type AbstractSystem end

# ==========================================================
# 1. Double Integrator Model
# ==========================================================
"""
    DoubleIntegrator

Represents a double-integrator system `x'' = u + d`. The safe region is a rectangle.
"""
struct DoubleIntegrator <: AbstractSystem
    dt::Float64
    safe_region::NamedTuple
end

function is_safe(system::DoubleIntegrator, state::Vector{Float64})
    x, v = state[1], state[2]
    xr = system.safe_region.x_lims
    vr = system.safe_region.v_lims
    return (xr[1] <= x <= xr[2]) && (vr[1] <= v <= vr[2])
end

function apply_dynamics(system::DoubleIntegrator, state::Vector{Float64}, action::Float64, disturbance::Float64)
    if !is_safe(system, state)
        return state
    end
    x, v = state[1], state[2]
    x_next = x + v * system.dt
    v_next = v + (action + disturbance) * system.dt
    return [x_next, v_next]
end


# ==========================================================
# 2. Inverted Pendulum Model
# ==========================================================
"""
    InvertedPendulum

Represents an inverted pendulum. The safe region is a ring between two ellipses.
"""
struct InvertedPendulum <: AbstractSystem
    dt::Float64
    g::Float64 # Gravity
    L::Float64 # Length of pendulum
    safe_region::NamedTuple # (a_outer, b_outer, a_inner, b_inner, center)
end

function is_safe(system::InvertedPendulum, state::Vector{Float64})
    theta, theta_dot = state[1], state[2]
    xc, yc = system.safe_region.center
    a_out, b_out = system.safe_region.a_outer, system.safe_region.b_outer
    a_in, b_in = system.safe_region.a_inner, system.safe_region.b_inner
    
    # Check if point is inside the larger ellipse
    inside_outer = ((theta - xc)^2 / a_out^2) + ((theta_dot - yc)^2 / b_out^2) <= 1
    # Check if point is inside the smaller ellipse (the hole)
    inside_inner = ((theta - xc)^2 / a_in^2) + ((theta_dot - yc)^2 / b_in^2) <= 1
    
    return inside_outer && !inside_inner
end

function apply_dynamics(system::InvertedPendulum, state::Vector{Float64}, action::Float64, disturbance::Float64)
    if !is_safe(system, state)
        return state
    end
    theta, theta_dot = state[1], state[2]
    
    theta_ddot = (system.g / system.L) * sin(theta) + (action + disturbance)
    
    theta_next = theta + theta_dot * system.dt
    theta_dot_next = theta_dot + theta_ddot * system.dt
    
    return [theta_next, theta_dot_next]
end


# ==========================================================
# 3. Dubins Car Model
# ==========================================================
"""
    DubinsCar

Represents a Dubins car model with constant forward velocity.
"""
struct DubinsCar <: AbstractSystem
    dt::Float64
    V::Float64 # Constant forward speed
    # For this model, let's assume the safe region is the same elliptical ring.
    safe_region::NamedTuple # (a_outer, b_outer, a_inner, b_inner, center)
end

function is_safe(system::DubinsCar, state::Vector{Float64})
    # NOTE: The safe region for Dubins car is usually on x and y, not theta.
    # We will check the safety on the first two state variables (x, y).
    x, y = state[1], state[2]
    xc, yc = system.safe_region.center
    a_out, b_out = system.safe_region.a_outer, system.safe_region.b_outer
    a_in, b_in = system.safe_region.a_inner, system.safe_region.b_inner

    inside_outer = ((x - xc)^2 / a_out^2) + ((y - yc)^2 / b_out^2) <= 1
    inside_inner = ((x - xc)^2 / a_in^2) + ((y - yc)^2 / b_in^2) <= 1

    return inside_outer && !inside_inner
end

function apply_dynamics(system::DubinsCar, state::Vector{Float64}, action::Float64, disturbance::Float64)
    if !is_safe(system, state)
        return state
    end
    x, y, theta = state[1], state[2], state[3]
    
    x_next = x + system.V * cos(theta) * system.dt
    y_next = y + system.V * sin(theta) * system.dt # Original Dubins car uses sin(theta) here
    theta_next = theta + (action + disturbance) * system.dt
    
    return [x_next, y_next, theta_next]
end


end # end of module Systems