module Plotting

using Plots
using LaTeXStrings
using ..Systems
using ..StateSpaces

gr() # Set the GR backend

export generate_and_save_plot

function generate_and_save_plot(;
    filepath::String,
    system::AbstractSystem,
    space::StateSpace,
    x1_params::Tuple,
    x2_params::Tuple,
    disturbance_props::NamedTuple,
    avr_results::Union{NamedTuple, Nothing} = nothing,
    mdr_results::Union{NamedTuple, Nothing} = nothing
)
    println("--- Generating final plot with manual annotation for title ---")
    
    # --- 1. Extract Data ---
    x1_coords = LinRange(x1_params[1], x1_params[2], x1_params[3])
    x2_coords = LinRange(x2_params[1], x2_params[2], x2_params[3])
    
    # --- 2. Initialize Plot (with NO title) ---
    p = plot(
        xlabel=L"x_1",
        ylabel=L"x_2",
        legend=:topright,
        framestyle=:box,
        tick_direction=:in,
        colorbar=false,
    )

    # --- 3. Plot Data (with empty labels) ---
    k1_min, k1_max = system.safe_region.x_lims
    k2_min, k2_max = system.safe_region.v_lims
    shape_x = [k1_min, k1_max, k1_max, k1_min, k1_min]
    shape_y = [k2_min, k2_min, k2_max, k2_max, k2_min]
    plot!(p, shape_x, shape_y,
        color=:black, linewidth=2.5, linestyle=:dash, label=""
    )

    if !isnothing(avr_results)
        contour!(p, x1_coords, x2_coords, avr_results.g_map', levels=[1.0], color=:blue, linewidth=2, label="")
    end
    if !isnothing(mdr_results)
        contour!(p, x1_coords, x2_coords, mdr_results.Z_map', levels=[0.0], color=:red, linewidth=2, label="")
    end

    # --- 4. Manually Add All Legend Entries ---
    plot!(p, [], [], color=:black, linewidth=2.5, linestyle=:dash, label="Constraint Set")
    if !isnothing(mdr_results)
        plot!(p, [], [], color=:red, linewidth=2, label=L"\{x \mid Z(x) = 0\}")
    end
    if !isnothing(avr_results)
        plot!(p, [], [], color=:blue, linewidth=2, label=L"\{x \mid g(x) = 1\}")
    end

    # --- 5. Manually Add the Title Using Annotate ---
    mu = disturbance_props.params.mean
    sigma = disturbance_props.params.sigma
    dmin, dmax = disturbance_props.bounds
    title_line1 = LaTeXString("\\mathcal{D} \\sim \\text{Normal}(\\mu=$(mu), \\sigma=$(sigma))")
    title_line2 = LaTeXString("d \\in [$(dmin), $(dmax)]")
    final_title = title_line1 * "\n" * title_line2

    # Calculate position for the title (top-center)
    title_x = (x1_coords[1] + x1_coords[end]) / 2
    # THIS IS THE CORRECTED LINE:
    title_y = Plots.ylims(p)[2] + 0.1 * (Plots.ylims(p)[2] - Plots.ylims(p)[1])

    annotate!(p, title_x, title_y, text(final_title, 12, :center))

    # --- 6. Save the Plot ---
    savefig(p, filepath)
    println("Plot saved to: $(filepath)")
end

end # end module