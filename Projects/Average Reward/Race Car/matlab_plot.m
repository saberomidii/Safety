clear; clc; close all;

% LaTeX styling
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultAxesTickLabelInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0, 'DefaultColorbarTickLabelInterpreter', 'latex');

% --- Grid Setup ---
x_1_min = -0.5;  x_1_max =  0.5;
x_2_min = -1.0;  x_2_max =  1.0;
num_points_state_1 = 62;
num_points_state_2 = 62;
num_points_state_3 = 16;

% --- Ellipse Parameters (updated) ---
a_outer = 0.4;  % Outer ellipse semi-major axis (x)
b_outer = 0.9;  % Outer ellipse semi-minor axis (y)

a_inner = 0.2;  % Inner ellipse semi-major axis (x)
b_inner = 0.7;  % Inner ellipse semi-minor axis (y)

% --- Load Data ---
G_Average = readmatrix('results_primal_1.0/G.csv');
G_Average = reshape(G_Average, [num_points_state_1, num_points_state_2, num_points_state_3]);

% --- Grid Mesh ---
x1_avg = linspace(x_1_min, x_1_max, num_points_state_1);
x2_avg = linspace(x_2_min, x_2_max, num_points_state_2);
[X1, X2] = meshgrid(x1_avg, x2_avg);

% --- Safe Set: single elliptical ring centered at (0,0) ---
inside_outer = ((X1).^2)/(a_outer^2) + ((X2).^2)/(b_outer^2) <= 1;
inside_inner = ((X1).^2)/(a_inner^2) + ((X2).^2)/(b_inner^2) <= 1;
in_safe_set = inside_outer & ~inside_inner;
mask = double(in_safe_set);

% Apply mask
G_Average_masked = G_Average .* mask;
G_Average_masked(mask == 0) = NaN;

% --- Colormap and Green for Safe Region ---
dark_green = [0, 0.5, 0];
cmap_length = 256;
custom_cmap = zeros(cmap_length, 3);
for i = 1:cmap_length
    if i <= cmap_length/2
        custom_cmap(i, :) = [1, (i-1)/(cmap_length/2-1), 0];      % Red to yellow
    else
        custom_cmap(i, :) = [1 - (i-cmap_length/2)/(cmap_length/2), 1, 0]; % Yellow to green
    end
end

% Prepare theta values
theta_list = linspace(-pi, pi, 16);

figure;
set(gcf, 'Position', [100, 100, 1200, 1000]);

for index_theta = 1:16
    % Extract the current slice for this theta
    G_slice = G_Average_masked(:, :, index_theta);
    
    % Create subplot 4x4 grid
    subplot(4,4,index_theta);
    hold on;
    
    % Plot contour
    colormap(custom_cmap);
    contour(x1_avg, x2_avg, G_slice, 100, 'LineWidth', 1.2);
    caxis([0 0.9]);
    
    % Colorbar only on last subplot for clarity
    if index_theta == 16
        cb = colorbar;
        ylabel(cb, 'Safety', 'FontSize', 14, 'Interpreter', 'none');
    end
    
    % Fill Safe Set
    high_safety_Avg = G_slice >= 0.99;
    [c_Avg, ~] = contour(x1_avg, x2_avg, high_safety_Avg, [0.5 0.5], '-k', 'LineWidth', 1.5);
    if ~isempty(c_Avg)
        contour_data = c_Avg;
        idx = 1;
        while idx < size(contour_data, 2)
            num_pts = contour_data(2, idx);
            fill(contour_data(1, idx+1:idx+num_pts), ...
                 contour_data(2, idx+1:idx+num_pts), ...
                 dark_green, 'EdgeColor', 'none');
            idx = idx + num_pts + 1;
        end
    end
    
    % Plot Ellipses (2D)
    theta_ellipse = linspace(0, 2*pi, 100);
    xe_outer = a_outer * cos(theta_ellipse);
    ye_outer = b_outer * sin(theta_ellipse);
    xe_inner = a_inner * cos(theta_ellipse);
    ye_inner = b_inner * sin(theta_ellipse);
    
    plot(xe_outer, ye_outer, 'k--', 'LineWidth', 2);
    plot(xe_inner, ye_inner, 'r--', 'LineWidth', 1.5);
    
    % Labels and Legend (only on first subplot)
    if index_theta == 1
        patch_handle = patch(NaN, NaN, dark_green);
        legend(patch_handle, {'Safe Set'}, 'Location', 'Northeast', ...
               'Interpreter', 'none', 'FontSize', 12);
    end
    
    % Axis labels with plain text
    xlabel('x_1', 'Interpreter', 'none', 'FontSize', 14);
    ylabel('x_2', 'Interpreter', 'none', 'FontSize', 14);
    
    % Title with theta in radians and degrees, no LaTeX
    theta_val = theta_list(index_theta);
    theta_deg = rad2deg(theta_val);
    title(sprintf('theta = %.2f rad, %.1f deg', theta_val, theta_deg), ...
          'Interpreter', 'none', 'FontSize', 14);
    
    % Axes Formatting
    x_ticks = linspace(x_1_min, x_1_max, 5);
    y_ticks = linspace(x_2_min, x_2_max, 5);
    set(gca, 'XTick', x_ticks, 'YTick', y_ticks);
    set(gca, 'XTickLabel', arrayfun(@(x) num2str(x, '%.1f'), x_ticks, 'UniformOutput', false));
    set(gca, 'YTickLabel', arrayfun(@(y) num2str(y, '%.1f'), y_ticks, 'UniformOutput', false));
    axis([x_1_min x_1_max x_2_min x_2_max]);
    axis square;
    grid on;
    hold off;
end
