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
        Z_lambda0_0 = readmatrix('MDR_Z_map_lambda_0.0.csv');
        Z_lambda0_01 = readmatrix("MDR_Z_map_lambda_0.01.csv");
        Z_lambda0_1 = readmatrix("MDR_Z_map_lambda_0.1.csv");
        Z_lambda0_2 = readmatrix("MDR_Z_map_lambda_0.2.csv");
        G_average = readmatrix("AVR_gain_map.csv");
    catch ME
        error('Failed to load data files. Make sure the CSV files are in the correct path.\n%s', ME.message);
    end
    
    % ---------------------------
    % 2. Reconstruct Grid & Define Parameters
    % ---------------------------
    % X is Position (Cols), V is Velocity (Rows) in the original matrix structure
    num_points_position = size(G_average, 2);
    num_points_velocity = size(G_average, 1);
    
    position_coords = linspace(-0.5, 0.5, num_points_position);
    velocity_coords = linspace(-1, 1, num_points_velocity);
    
    constraint_position_lim = [-0.3, 0.3];
    constraint_velocity_lim = [-0.6, 0.6];
    
    % ###########################################################################
    % ## FIGURE: Safe Set Comparison
    % ###########################################################################
    
    fig = figure('Position', [100, 100, 700, 650]);
    ax = axes('Parent', fig);
    hold(ax, 'on');
    
    % --- Plot Data ---
    % CRITICAL: Pass limits as [X_lim, Y_lim] -> [Velocity_lim, Position_lim]
    plotConstraintSet(ax, constraint_velocity_lim, constraint_position_lim);
    
    % CRITICAL: Transpose data matrices (.') and swap coordinate vectors (velocity, position)
    contour(ax, velocity_coords, position_coords, G_average.', [1.0, 1.0], 'LineWidth', 3.5, 'Color', 'k', 'LineStyle', '-', 'DisplayName', 'AVR Safe Set ($g^*(s)=1$)');
    contour(ax, velocity_coords, position_coords, Z_lambda0_0.', [0.0 0.0], 'LineWidth', 2.5, 'Color', [0.8500 0.3250 0.0980], 'LineStyle', '-', 'DisplayName', 'MDR ($Z(x)=0, \lambda=0.0$)');
    
    % Data masking is kept on the original matrix structure (Rows=V, Cols=X). 
    % When transposed for plotting, the masked area remains visually correct.
    Z_lambda0_01(51:135,1:55) = -10;
    Z_lambda0_01(65:151,148:201) = -10;
    contour(ax, velocity_coords, position_coords, Z_lambda0_01.', [0.0, 0.0], 'LineWidth', 2, 'Color', [0 0.4470 0.7410], 'LineStyle', '--', 'DisplayName', 'MDR ($Z(x)=0, \lambda=0.1$)');
    
    Z_lambda0_1(51:135,1:55) = -10;
    Z_lambda0_1(65:151,148:201) = -10;
    contour(ax, velocity_coords, position_coords, Z_lambda0_1.', [0.0, 0.0], 'LineWidth', 2, 'Color', [0.4660 0.6740 0.1880], 'LineStyle', ':', 'DisplayName', 'MDR ($Z(x)=0, \lambda=0.2$)');
    
    Z_lambda0_2(51:135,1:55) = -10;
    Z_lambda0_2(65:151,148:201) = -10;
    contour(ax, velocity_coords, position_coords, Z_lambda0_2.', [0.0,0.0], 'LineWidth', 2, 'Color', [0.4940 0.1840 0.5560], 'LineStyle', '-.', 'DisplayName', 'MDR ($Z(x)=0, \lambda=0.3$)');
    
    hold(ax, 'off');
    
    % --- Format the plot for the paper ---
    formatAxesAndLabels(ax);
    
    % ---------------------------
    % 3. Save Figure to PDF without extra whitespace
    % ---------------------------
    output_filename = 'Safe_Set_Comparison_Rotated_Set2.pdf';
    fprintf('Saving figure to PDF...\n');
    
    set(ax, 'LooseInset', get(ax, 'TightInset'));
    
    exportgraphics(fig, output_filename, 'ContentType', 'vector', 'BackgroundColor', 'white');
    fprintf('Figure saved as %s\n', output_filename);
end
% ###########################################################################
% ## Helper Functions
% ###########################################################################
function plotConstraintSet(ax, velocity_lim, position_lim)
    % Plots a dashed box representing the constraint set C on the given axes.
    % The input limits are now interpreted as: X-axis (Velocity), Y-axis (Position)
    
    box_x = [velocity_lim(1), velocity_lim(2), velocity_lim(2), velocity_lim(1), velocity_lim(1)]; % Velocity corners
    box_y = [position_lim(1), position_lim(1), position_lim(2), position_lim(2), position_lim(1)]; % Position corners
    
    % Plotting (Velocity, Position)
    plot(ax, box_x, box_y, 'Color', [0.5 0.5 0.5], 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', 'Constraint Set ($\mathcal{C}$)');
end
function formatAxesAndLabels(ax)
    % Sets up the title, labels, grid, and legend with paper-consistent notation.
    grid(ax, 'on');
    axis(ax, 'equal');
    
    % CRITICAL: Swap axis limits for the rotation
    xlim(ax, [-0.65, 0.65]); % New X-axis (Velocity) limits
    ylim(ax, [-0.4, 0.4]);    % New Y-axis (Position) limits
    
    % Set font size for the tick labels on both axes
    set(ax, 'FontSize', 16);
    title(ax, '', 'Interpreter', 'latex', 'FontSize', 18);
    
    % CRITICAL: Swap labels
    xlabel(ax, 'Velocity ($\mathrm{v}$)', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel(ax, 'Position ($\mathrm{x}$)', 'Interpreter', 'latex', 'FontSize', 16);
    
    lgd = legend(ax, 'show');
    
    % --- MANUAL POSITIONING FOR LEGEND ---
    % 1. Get the current position vector [left, bottom, width, height]
    legend_pos = get(lgd, 'Position');
    
    % 2. Define the new (left, bottom) coordinates (normalized to figure size, 0 to 1)
    new_left = 0.1;  % Move it slightly right (5% from left edge of figure)
    new_bottom = 0.275; % Move it up (30% from bottom edge of figure)
    
    % 3. Set the new position using the calculated width/height
    set(lgd, 'Position', [new_left, new_bottom, legend_pos(3), legend_pos(4)], ...
        'Interpreter', 'latex', 'FontSize', 14); 
    % Note: If you use the 'Position' property, the 'Location' property is ignored/overwritten.
end