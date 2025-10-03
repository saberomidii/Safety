module CostFunctions

export make_box_cost_function, make_ellipse_cost_function

"""
    make_box_cost_function(safe_region)

Returns a function that calculates the signed distance to a rectangular box.
"""
function make_box_cost_function(safe_region::NamedTuple)
    
    # This is the actual cost function that will be returned and used
    function cost_function(state::Vector{Float64})
        px, py = state[1], state[2]
        k1_min, k1_max = safe_region.x_lims
        k2_min, k2_max = safe_region.v_lims

        dx = max(k1_min - px, 0.0, px - k1_max)
        dy = max(k2_min - py, 0.0, py - k2_max)
        
        dist_outside = sqrt(dx^2 + dy^2)
        return dist_outside > 0 ? dist_outside : -min(px - k1_min, k1_max - px, py - k2_min, k2_max - py)
    end
    
    return cost_function
end

"""
    make_ellipse_cost_function(safe_region)

Returns a function that calculates an approximate signed distance to an elliptical ring.
"""
function make_ellipse_cost_function(safe_region::NamedTuple)

    function cost_function(state::Vector{Float64})
        x, y = state[1], state[2]
        xc, yc = safe_region.center
        a_out, b_out = safe_region.a_outer, safe_region.b_outer
        a_in, b_in = safe_region.a_inner, safe_region.b_inner

        val_outer = ((x - xc)^2 / a_out^2) + ((y - yc)^2 / b_out^2)
        val_inner = ((x - xc)^2 / a_in^2) + ((y - yc)^2 / b_in^2)

        if val_outer > 1.0
            return sqrt(val_outer) - 1.0
        elseif val_inner < 1.0
            return 1.0 - sqrt(val_inner)
        else
            return -min(sqrt(val_outer) - 1.0, 1.0 - sqrt(val_inner))
        end
    end
    
    return cost_function
end

end # end module