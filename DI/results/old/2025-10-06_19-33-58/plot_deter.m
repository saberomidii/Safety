function plot_deter()
    clear; clc;
    close all;
    fprintf('Plotting results from Julia simulation...\n');
    
    % ---------------------------
    % 1. Load Data
    % ---------------------------
    % Data is transposed here upon loading
    Z_0 = readmatrix('MDR_Z_map_lambda_0.0.csv');
    Z_1 = readmatrix("MDR_Z_map_lambda_0.01.csv");
    Z_2 = readmatrix("MDR_Z_map_lambda_0.1.csv");
    Z_3 = readmatrix("MDR_Z_map_lambda_0.2.csv");
    G_average  = readmatrix("AVR_gain_map.csv");
        
    % ---------------------------
    % 2. Reconstruct Grid & Parameters
    % ---------------------------
    NUM_POINTS_STATE_1 = size(G_average, 1);
    NUM_POINTS_STATE_2 = size(G_average, 2);
    x1 = linspace(-1, 5, NUM_POINTS_STATE_1);
    x2 = linspace(-5, 5, NUM_POINTS_STATE_2);
    
    c_min_1 = 0.0; c_max_1 = 4.0;
    c_min_2 = -3.0; c_max_2 = 3.0;
    
    % ###########################################################################
    % ## FIGURE: RESULTS VISUALIZATION
    % ###########################################################################
    
    fig = figure('Position', [100, 100, 800, 700]); 
    hold on;
    
    plotConstraintSetBox(c_min_1, c_max_1, c_min_2, c_max_2);
    
    % --- CORRECTION: Removed extra transpose (') from plotting calls ---
    contour(x1, x2, G_average, [1.0, 1.0], "LineWidth", 3, "Color", 'k', 'LineStyle', '-', 'DisplayName', 'AVR Safe Set (g=1)');
    contour(x1, x2, Z_0, [0.0 0.0], "LineWidth", 2.5, "Color", 'r', 'LineStyle', '-', 'DisplayName', 'MDR ($\lambda=0.0$)');
    contour(x1, x2, Z_1, [0.0, 0.0], "LineWidth", 2, "Color", 'b', 'LineStyle', '--', 'DisplayName', 'MDR ($\lambda=0.01$)');
    contour(x1, x2, Z_2, [0.0, 0.0], "LineWidth", 2, "Color", 'g', 'LineStyle', ':', 'DisplayName', 'MDR ($\lambda=0.1$)');
    contour(x1, x2, Z_3, [0.0, 0.0], "LineWidth", 2, "Color", 'm', 'LineStyle', '-.', 'DisplayName', 'MDR ($\lambda=0.2$)');
    
    hold off;
    
    setupAxes();
    
    fprintf('\nSaving figure to PDF...\n');
    exportgraphics(fig, 'AVR_vs_MDR_plot.pdf', 'ContentType', 'vector', 'BackgroundColor', 'white');
    fprintf('Figure saved as AVR_vs_MDR_plot.pdf\n');

end

% ---------------------------
% Helper Functions
% ---------------------------
function plotConstraintSetBox(x_min, x_max, y_min, y_max)
    box_x = [x_min, x_max, x_max, x_min, x_min];
    box_y = [y_min, y_min, y_max, y_max, y_min];
    plot(box_x, box_y, 'Color', [0.4 0.4 0.4], 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', 'Constraint Set');
end

function setupAxes()
    grid on; 
    axis equal;
    
    xlim([-1.5, 5.5]);
    ylim([-5.5, 5.5]);
    
    % --- UPDATED: New title with Normal Distribution information ---
    title_text = 'AVR vs. MDR for $\mathcal{D} \sim \mathcal{N}(\mu=0, \sigma=1)$, $d \in [-1,1]$';
    
    title(title_text, 'Interpreter', 'latex', 'FontSize', 16);
    xlabel('Position ($x_1$)', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('Velocity ($x_2$)', 'Interpreter', 'latex', 'FontSize', 14);
    
    lgd = legend('show');
    set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'FontSize', 12);
end