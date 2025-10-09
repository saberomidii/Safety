function generate_paper_figure()
    % This function loads simulation results, generates a comparative plot of the 
    % Average Reward (AVR) and Minimum Discounted Reward (MDR) safe sets,
    % and saves the final figure as a high-quality, tightly cropped PDF.

    clear; clc; close all;
    fprintf('Generating the safe set comparison figure for the paper...\n');
    
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
    num_points_x = size(G_average, 2); 
    num_points_v = size(G_average, 1); 
    
    x_coords = linspace(-0.5, 0.5, num_points_x);
    v_coords = linspace(-1.0,1.0, num_points_v);
    
    constraint_x_lim = [-0.3, 0.3];
    constraint_v_lim = [-0.6, 0.6];
    
    % ###########################################################################
    % ## FIGURE: Safe Set Comparison
    % ###########################################################################
    
    fig = figure('Position', [100, 100, 700, 650]); 
    ax = axes('Parent', fig);
    hold(ax, 'on');
    
    % --- Plot Data ---
    plotConstraintSet(ax, constraint_x_lim, constraint_v_lim);
    epsilon =0.000;
    contour(ax, x_coords, v_coords, G_average, [1, 1], 'LineWidth', 3.5, 'Color', 'k', 'LineStyle', '-', 'DisplayName', 'AVR Safe Set ($g^*(s)=1$)');
    contour(ax, x_coords, v_coords, Z_lambda0_0, [epsilon epsilon], 'LineWidth', 2.5, 'Color', [0.8500 0.3250 0.0980], 'LineStyle', '-', 'DisplayName', 'MDR ($Z(x)=0, \lambda=0.0$)');
    contour(ax, x_coords, v_coords, Z_lambda0_01, [epsilon, epsilon], 'LineWidth', 2, 'Color', [0 0.4470 0.7410], 'LineStyle', '--', 'DisplayName', 'MDR ($Z(x)=0, \lambda=0.05$)');
    contour(ax, x_coords, v_coords, Z_lambda0_1, [epsilon, epsilon], 'LineWidth', 2, 'Color', [0.4660 0.6740 0.1880], 'LineStyle', ':', 'DisplayName', 'MDR ($Z(x)=0, \lambda=0.1$)');
    contour(ax, x_coords, v_coords, Z_lambda0_2, [epsilon, epsilon], 'LineWidth', 2, 'Color', [0.4940 0.1840 0.5560], 'LineStyle', '-.', 'DisplayName', 'MDR ($Z(x)=0, \lambda=0.15$)');
    
    hold(ax, 'off');
    
    % --- Format the plot for the paper ---
    formatAxesAndLabels(ax);
    
    % ---------------------------
    % 3. Save Figure to PDF without extra whitespace
    % ---------------------------
    output_filename = 'Safe_Set_Comparison.pdf';
    fprintf('Saving figure to PDF...\n');
    
    set(ax, 'LooseInset', get(ax, 'TightInset'));
    
    exportgraphics(fig, output_filename, 'ContentType', 'vector', 'BackgroundColor', 'white');
    fprintf('Figure saved as %s\n', output_filename);
end

% ###########################################################################
% ## Helper Functions
% ###########################################################################

function plotConstraintSet(ax, x_lim, v_lim)
    % Plots a dashed box representing the constraint set C on the given axes.
    box_x = [x_lim(1), x_lim(2), x_lim(2), x_lim(1), x_lim(1)];
    box_v = [v_lim(1), v_lim(1), v_lim(2), v_lim(2), v_lim(1)];
    plot(ax, box_x, box_v, 'Color', [0.5 0.5 0.5], 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', 'Constraint Set ($\mathcal{C}$)');
end

function formatAxesAndLabels(ax)
    % Sets up the title, labels, grid, and legend with paper-consistent notation.
    grid(ax, 'on'); 
    axis(ax, 'equal');
    
    % UPDATED: Axis limits have been adjusted to "zoom in".
    xlim(ax, [-0.5, 0.5]);
    ylim(ax, [-1, 0.65]);
    
    % Set font size for the tick labels on both axes
    set(ax, 'FontSize', 16);

    title(ax, '', 'Interpreter', 'latex', 'FontSize', 18);
    
    xlabel(ax, 'Position ($\mathrm{x}$)', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel(ax, 'Velocity ($\mathrm{v}$)', 'Interpreter', 'latex', 'FontSize', 16);
    
    lgd = legend(ax, 'show');
    set(lgd, 'Interpreter', 'latex', 'Location', 'southwest', 'FontSize', 14);
end