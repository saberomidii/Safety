function generate_paper_figure()
    % This function loads simulation results, generates a comparative plot of the 
    % Average Reward (AVR) and Minimum Discounted Reward (MDR) safe sets,
    % with VELOCITY on the X-axis and POSITION on the Y-axis,
    % and saves the final figure as a high-quality, tightly cropped PDF.
    clear; clc; close all;
    fprintf('Generating the safe set comparison figure for the paper (Rotated Axes)...\n');
    
    % ---------------------------
    % 1. Load Simulation Data
    % ---------------------------
    try
        Z_lambda0_0  = readmatrix('MDR_Z_map_lambda_0.0.csv');
        Z_lambda0_01 = readmatrix("MDR_Z_map_lambda_0.01.csv");
        Z_lambda0_1  = readmatrix("MDR_Z_map_lambda_0.1.csv");
        Z_lambda0_2  = readmatrix("MDR_Z_map_lambda_0.2.csv");
        G_average    = readmatrix("AVR_gain_map.csv");
    catch ME
        error('Failed to load data files. Make sure the CSV files are in the correct path.\n%s', ME.message);
    end
        
    % ---------------------------
    % 2. Reconstruct Grid & Define Parameters
    % ---------------------------
    % Position (X-dim/columns) and Velocity (Y-dim/rows) from the input data
    num_points_position = size(G_average, 2); 
    num_points_velocity = size(G_average, 1); 
    
    position_coords = linspace(-1, 5, num_points_position);
    velocity_coords = linspace(-5, 5, num_points_velocity);
    
    constraint_position_lim = [0.0, 4.0];
    constraint_velocity_lim = [-3.0, 3.0];
    
    % ###########################################################################
    % ## FIGURE: Safe Set Comparison
    % ###########################################################################
    
    fig = figure('Position', [100, 100, 700, 650]); 
    ax = axes('Parent', fig);
    hold(ax, 'on');
    
    % --- Plot Data ---
    % NOTE: The order of constraint limits is swapped in the call: [X_lim, Y_lim] -> [V_lim, X_lim]
    plotConstraintSet(ax, constraint_velocity_lim, constraint_position_lim);
    
    epsilon =0.015;
    
    % CRITICAL CHANGE: The coordinate vectors are swapped, and the data matrix is transposed.
    % contour(ax, X_coords, Y_coords, Data_TRANSPOSED, ...)
    contour(ax, velocity_coords, position_coords, G_average.', [1.0, 1.0], 'LineWidth', 3.5, 'Color', 'k', 'LineStyle', '-', 'DisplayName', 'AVR Safe Set ($g^*(s)=1$)');
    contour(ax, velocity_coords, position_coords, Z_lambda0_0.', [0.0 0.0], 'LineWidth', 3.5, 'Color', [0.8500 0.3250 0.0980], 'LineStyle', '-', 'DisplayName', 'MDR ($Z(x)=0, \lambda=0.0$)');
    contour(ax, velocity_coords, position_coords, Z_lambda0_01.', [epsilon, epsilon], 'LineWidth', 3, 'Color', [0 0.4470 0.7410], 'LineStyle', '--', 'DisplayName', 'MDR ($Z(x)=0, \lambda=0.01$)');
    contour(ax, velocity_coords, position_coords, Z_lambda0_1.', [epsilon, epsilon], 'LineWidth', 3, 'Color', [0.4660 0.6740 0.1880], 'LineStyle', ':', 'DisplayName', 'MDR ($Z(x)=0, \lambda=0.015$)');
    contour(ax, velocity_coords, position_coords, Z_lambda0_2.', [epsilon,epsilon], 'LineWidth', 3, 'Color', [0.4940 0.1840 0.5560], 'LineStyle', '-.', 'DisplayName', 'MDR ($Z(x)=0, \lambda=0.02$)');
    
    hold(ax, 'off');
    
    % --- Format the plot for the paper ---
    formatAxesAndLabels(ax);
    
    % ---------------------------
    % 3. Save Figure to PDF without extra whitespace
    % ---------------------------
    output_filename = 'Safe_Set_Comparison_Rotated.pdf';
    fprintf('Saving figure to PDF...\n');
    
    set(ax, 'LooseInset', get(ax, 'TightInset'));
    
    exportgraphics(fig, output_filename, 'ContentType', 'vector', 'BackgroundColor', 'white');
    fprintf('Figure saved as %s\n', output_filename);
end
% ###########################################################################
% ## Helper Functions
% ###########################################################################
function plotConstraintSet(ax, x_lim, y_lim)
    % Plots a dashed box representing the constraint set C on the given axes.
    % The input limits are now interpreted as: x_lim = Velocity, y_lim = Position
    box_x = [x_lim(1), x_lim(2), x_lim(2), x_lim(1), x_lim(1)]; % Velocity corners
    box_y = [y_lim(1), y_lim(1), y_lim(2), y_lim(2), y_lim(1)]; % Position corners
    
    % CRITICAL CHANGE: plotting (Velocity, Position)
    plot(ax, box_x, box_y, 'Color', [0.5 0.5 0.5], 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', 'Constraint Set ($\mathcal{C}$)');
end
function formatAxesAndLabels(ax)
    % Sets up the title, labels, grid, and legend with paper-consistent notation.
    grid(ax, 'on'); 
    axis(ax, 'equal');
    
    % CRITICAL CHANGE: Swapping xlim and ylim values and labels
    
    % Old X (Position) limits: [-0.5, 4.5] -> New Y (Position) limits
    % Old Y (Velocity) limits: [-4.75, 3.5] -> New X (Velocity) limits
    xlim(ax, [-4.75, 3.5]); % New X-axis (Velocity) limits
    ylim(ax, [-0.5, 4.5]);  % New Y-axis (Position) limits
    
    % Set font size for the tick labels on both axes
    set(ax, 'FontSize', 16);
    title(ax, '', 'Interpreter', 'latex', 'FontSize', 18);
    
    % CRITICAL CHANGE: Swapping labels
    xlabel(ax, 'Velocity ($\mathrm{v}$)', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel(ax, 'Position ($\mathrm{x}$)', 'Interpreter', 'latex', 'FontSize', 16);
    
    lgd = legend(ax, 'show');
    set(lgd, 'Interpreter', 'latex', 'Location', 'southwest', 'FontSize', 14);
end