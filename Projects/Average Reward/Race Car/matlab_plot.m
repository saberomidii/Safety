clear; clc; close all;

% LaTeX styling
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultAxesTickLabelInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0, 'DefaultColorbarTickLabelInterpreter', 'latex');

% --- Grid Setup ---
x_1_min = -0.5;  x_1_max =  0.5;
x_2_min = -1.0;  x_2_max =  1.0;
num_points_state_1 = 30;
num_points_state_2 = 30;
num_points_state_3 = 30;

% --- Ellipse Parameters (updated) ---
a_outer = 0.4;  % Outer ellipse semi-major axis (x)
b_outer = 0.9;  % Outer ellipse semi-minor axis (y)

a_inner = 0.2;  % Inner ellipse semi-major axis (x)
b_inner = 0.7;  % Inner ellipse semi-minor axis (y)

% --- Load Data ---
G_Average = readmatrix('results_primal_2/G.csv');
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

% === 2D CONTOUR PLOT ===
G_slice = G_Average_masked(:, :, 22);
figure;
set(gcf, 'Position', [100, 100, 600, 500]);
hold on;
colormap(custom_cmap);
contour(x1_avg, x2_avg, G_slice, 100, 'LineWidth', 1.2);
caxis([0 0.9]);

% Colorbar
cb = colorbar;
ylabel(cb, 'Safety', 'FontSize', 14, 'Interpreter', 'latex');

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
theta = linspace(0, 2*pi, 100);
xe_outer = a_outer * cos(theta);
ye_outer = b_outer * sin(theta);
xe_inner = a_inner * cos(theta);
ye_inner = b_inner * sin(theta);

plot(xe_outer, ye_outer, 'k--', 'LineWidth', 2);
plot(xe_inner, ye_inner, 'r--', 'LineWidth', 1.5);

% Labels and Legend
patch_handle = patch(NaN, NaN, dark_green);
legend(patch_handle, {'Safe Set'}, 'Location', 'Northeast', ...
       'Interpreter', 'latex', 'FontSize', 12);

xlabel('$x_1$', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('$x_2$', 'FontSize', 20, 'FontWeight', 'bold');
title('Average Reward Safety Set', 'FontSize', 18, 'Interpreter', 'latex');

% Axes Formatting
x_ticks = linspace(x_1_min, x_1_max, 5);
y_ticks = linspace(x_2_min, x_2_max, 5);
set(gca, 'XTick', x_ticks, 'YTick', y_ticks);
set(gca, 'XTickLabel', arrayfun(@(x) ['$' num2str(x, '%.1f') '$'], x_ticks, 'UniformOutput', false));
set(gca, 'YTickLabel', arrayfun(@(y) ['$' num2str(y, '%.1f') '$'], y_ticks, 'UniformOutput', false));
axis([x_1_min x_1_max x_2_min x_2_max]);
axis square;
grid on;

% --- 3D SURFACE PLOT ---
figure;
set(gcf, 'Position', [100, 100, 1400, 600]);
tiledlayout(5, 6, 'Padding', 'compact', 'TileSpacing', 'compact');

for theta_idx = 1:num_points_state_3
    nexttile;
    slice_data = G_Average_masked(:, :, theta_idx);
    surf(X1, X2, slice_data, 'EdgeColor', 'none');
    view(45, 30);
    zlim([0 1]);
    caxis([0 0.9]);
    colormap(custom_cmap);
    
    hold on;
    % Plot outer ellipse on surface plot as black dashed line at z=0
    plot3(a_outer*cos(theta), b_outer*sin(theta), zeros(size(theta)), 'k--', 'LineWidth', 1.5);
    
    % Optional: plot inner ellipse as well
    plot3(a_inner*cos(theta), b_inner*sin(theta), zeros(size(theta)), 'r--', 'LineWidth', 1);
    
    title(['$\theta\_$' num2str(theta_idx)], 'Interpreter', 'latex', 'FontSize', 12);
    xlabel('$x_1$', 'Interpreter', 'latex');
    ylabel('$x_2$', 'Interpreter', 'latex');
    zlabel('$G$', 'Interpreter', 'latex');
    set(gca, 'FontSize', 10, 'TickLabelInterpreter', 'latex');
end

% Colorbar for all plots
cb = colorbar('Position', [0.93 0.1 0.015 0.8]);
ylabel(cb, 'Safety', 'FontSize', 14, 'Interpreter', 'latex');


% === 2D CONTOUR PLOTS FOR ALL THETA SLICES IN ONE FIGURE ===
figure;
set(gcf, 'Position', [100, 100, 1400, 700]);
tiledlayout(2, ceil(num_points_state_3/2), 'Padding', 'compact', 'TileSpacing', 'compact');

for theta_idx = 1:num_points_state_3
    nexttile;
    G_slice = G_Average_masked(:, :, theta_idx);
    
    contour(x1_avg, x2_avg, G_slice, 100, 'LineWidth', 1.2);
    caxis([0 0.9]);
    colormap(custom_cmap);
    hold on;
    
    % Plot ellipses
    plot(xe_outer, ye_outer, 'k--', 'LineWidth', 2);
    plot(xe_inner, ye_inner, 'r--', 'LineWidth', 1.5);
    
    xlabel('$x_1$', 'Interpreter', 'latex');
    ylabel('$x_2$', 'Interpreter', 'latex');
    title(['$\theta\_$' num2str(theta_idx)], 'Interpreter', 'latex', 'FontSize', 14);
    axis([x_1_min x_1_max x_2_min x_2_max]);
    axis square;
    grid on;
    
    % Optional: remove tick labels on inner plots for cleaner look
    if mod(theta_idx-1, ceil(num_points_state_3/2)) ~= 0
        set(gca, 'YTickLabel', []);
    end
    if theta_idx <= num_points_state_3 - ceil(num_points_state_3/2)
        set(gca, 'XTickLabel', []);
    end
end

% Colorbar for all subplots (single colorbar for whole figure)
cb = colorbar('Position', [0.93 0.1 0.015 0.8]);
ylabel(cb, 'Safety', 'FontSize', 14, 'Interpreter', 'latex');

