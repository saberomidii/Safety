function plot_avr_heatmap()
    % This function loads AVR data for the double integrator, generates a 
    % detailed heatmap that specifically highlights the one-level set, 
    % overlays the constraint set, and saves the figure as a PDF.
    
    % Setup environment
    clear; clc; close all;
    fprintf('Generating the AVR value function heatmap for the Double Integrator...\n');
    
    % ---------------------------
    % 1. Load Simulation Data
    % ---------------------------
    try
        G_average = readmatrix("AVR_gain_map.csv");
        fprintf('AVR_gain_map.csv loaded successfully.\n');
    catch ME
        error('Failed to load data file. Make sure "AVR_gain_map.csv" is in the correct path.\n%s', ME.message);
    end
        
    % ---------------------------
    % 2. Reconstruct Grid & Define Parameters
    % ---------------------------
    num_points_v = size(G_average, 1); 
    num_points_x = size(G_average, 2); 
    
    x_coords = linspace(-1, 5, num_points_x);
    v_coords = linspace(-5, 5, num_points_v);
    
    constraint_x_lim = [0.0, 4.0];
    constraint_v_lim = [-3.0, 3.0];
    
    % ###########################################################################
    % ## FIGURE: AVR Value Function Heatmap for Double Integrator
    % ###########################################################################
    
    fig = figure('Position', [100, 100, 700, 600]); 
    ax = axes('Parent', fig);
    hold(ax, 'on');
    
    % --- Create and plot the heatmap ---
    imagesc(ax, x_coords, v_coords, G_average);
    set(ax, 'YDir', 'normal');
    
    % --- Plot the constraint set on top ---
    plotConstraintSet(ax, constraint_x_lim, constraint_v_lim);
    
    hold(ax, 'off');
    
    % --- Format the plot for the paper ---
    formatAxesAndLabels(ax);
    
    % ---------------------------
    % 3. Save Figure to PDF
    % ---------------------------
    output_filename = 'AVR_Heatmap_DoubleIntegrator_Final.pdf';
    fprintf('Saving figure to PDF...\n');
    
    set(ax, 'LooseInset', get(ax, 'TightInset'));
    exportgraphics(fig, output_filename, 'ContentType', 'vector', 'BackgroundColor', 'white');
    fprintf('Figure saved as %s\n', output_filename);
end

% ###########################################################################
% ## Helper Functions
% ###########################################################################

function plotConstraintSet(ax, x_lim, v_lim)
    % Plots a solid black box representing the constraint set C.
    box_x = [x_lim(1), x_lim(2), x_lim(2), x_lim(1), x_lim(1)];
    box_v = [v_lim(1), v_lim(1), v_lim(2), v_lim(2), v_lim(1)];
    plot(ax, box_x, box_v, 'Color', 'k', 'LineWidth', 2.5, 'LineStyle', '-');
end

function formatAxesAndLabels(ax)
    % Sets up the title, labels, and a special colormap to highlight the g=1 level set.
    
    % --- UPDATED: Create a precise Red-to-Green colormap ---
    % This maps g=0 to darkest red and g=1 to darkest green, with a rich
    % gradient in between for values 0 < g < 1.
    
    % Part 1: Gradient for values between 0 and 1 (255 colors)
    value_points_gradient = [0.0, 0.2, 0.5, 0.8, 0.99];
    color_anchors_gradient = [
        0.5, 0.0, 0.0;  % Darkest Red (at g=0)
        1.0, 0.0, 0.0;  % Bright Red
        1.0, 0.8, 0.0;  % Yellow/Orange
        0.5, 1.0, 0.5;  % Light Green
        0.0, 0.8, 0.0   % Bright Green (just before g=1)
    ];
    gradient_map = interp1(value_points_gradient, color_anchors_gradient, linspace(0, 1, 255));
    
    % Part 2: Special color for the one-level set (g=1)
    darkest_green = [0.0, 0.4, 0.0];
    
    % Part 3: Combine them into the final colormap
    final_colormap = [gradient_map; darkest_green];
    colormap(ax, final_colormap);
    
    % Crucially, set the color axis limits to [0, 1]. This forces all values
    % of 1.0 to map to the very last color in our map (darkest_green).
    caxis(ax, [0 1]);
    
    % --- Add and format the color bar ---
    h_bar = colorbar(ax);
    ylabel(h_bar, 'Value ($g(s)$)', 'Interpreter', 'latex', 'FontSize', 16);
    
    % --- Set titles and labels ---
    title(ax, 'AVR Value Function (Double Integrator)', 'Interpreter', 'latex', 'FontSize', 18);
    xlabel(ax, 'Position ($\mathrm{x}$)', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel(ax, 'Velocity ($\mathrm{v}$)', 'Interpreter', 'latex', 'FontSize', 16);
    
    % --- Final visual adjustments ---
    grid(ax, 'off'); 
    axis(ax, 'tight');
    box(ax, 'on');
    set(ax, 'FontSize', 16, 'Layer', 'top'); 
end