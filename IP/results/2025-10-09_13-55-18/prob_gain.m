% MATLAB script to generate a single plot for G
clear all;
clc;
% close all;

% --- 1. Setup and Configuration ---

% Set LaTeX as the default interpreter for a professional look
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultAxesTickLabelInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0, 'DefaultColorbarTickLabelInterpreter', 'latex');

% Define font sizes for the plot
axis_font_size = 18;  % Font size for the tick labels
label_font_size = 20; % Font size for the x/y labels

% Define grid settings for the plot axes
x_1_min = -0.5;
x_1_max = 0.5;
x_2_min = -1.0;
x_2_max = 1.0;

% Define the rectangular safe set boundary for visualization
safe_x = [-0.3, 0.3, 0.3, -0.3, -0.3];
safe_v = [-0.6, -0.6, 0.6, 0.6, -0.6];
% --- 2. Data Loading and Preparation ---

% Load data from the specified CSV file
G = readmatrix('AVR_gain_map.csv');

% Create coordinate grid vectors corresponding to the dimensions of G
x1_grid = linspace(x_1_min, x_1_max, size(G, 2));
x2_grid = linspace(x_2_min, x_2_max, size(G, 1));

% Create a mask to highlight the G=0 boundary area
% This prevents plotting intermediate contours in the outer regions
g_mask = G;
g_mask(1:39, 1:201) = -1;
g_mask(163:201, 1:201) = -1;
g_mask(1:201, 1:39) = -1;
g_mask(1:201, 163:201) = -1;
% --- 3. Create the Custom Colormap with Orange and Yellow --- %%% UPDATED SECTION %%%
% Define the number of levels/colors needed for the gradient
num_levels = 9999;

% Define the key colors (nodes) for the gradient: Red -> Yellow -> Green
color_red    = [1, 0.65 0];
color_yellow = [1, 1, 0];
color_green  = [0.5, 1.0, 0.0];

% Calculate how many steps are in each segment
num_steps_part1 = round(num_levels / 2); % Steps for Red -> Yellow
num_steps_part2 = num_levels - num_steps_part1; % Steps for Yellow -> Green

% Create the first segment (Red to Yellow)
r1 = linspace(color_red(1), color_yellow(1), num_steps_part1);
g1 = linspace(color_red(2), color_yellow(2), num_steps_part1);
b1 = linspace(color_red(3), color_yellow(3), num_steps_part1);
part1 = [r1', g1', b1'];

% Create the second segment (Yellow to Green)
% We create N+1 points and then remove the first to avoid duplicating yellow
r2 = linspace(color_yellow(1), color_green(1), num_steps_part2 + 1);
g2 = linspace(color_yellow(2), color_green(2), num_steps_part2 + 1);
b2 = linspace(color_yellow(3), color_green(3), num_steps_part2 + 1);
part2 = [r2(2:end)', g2(2:end)', b2(2:end)'];

% Combine the segments into the final colormap
custom_cmap = [part1; part2];

% --- 4. Generate the Plot ---

% Create and configure the figure window
figure;
set(gcf, 'Position', [100, 100, 750, 550], 'Color', 'white');
hold on; % Keep all subsequent plots on the same axes

% Plot the specific filled contours for G=0 (unsafe) and G=1 (target) sets
contourf(x1_grid, x2_grid, g_mask, [0, 0], 'FaceColor', [1 0 0], 'EdgeColor', 'k'); % Grey for unsafe

contourf(x1_grid, x2_grid, g_mask, [0.9999, 0.9999], 'FaceColor', 'g', 'EdgeColor', 'k'); % Green for target

% Plot the 999 intermediate contour levels with the color gradient
for i = 1:num_levels
    level = 0.0001 * i;
    line_color = custom_cmap(i, :); % Select the i-th color
    contour(x1_grid, x2_grid, G, [level, level], 'LineWidth', 2.5, 'Color', line_color);
end

% Plot the rectangular boundary over the top
plot(safe_x, safe_v, 'k--', 'LineWidth', 3);


contourf(x1_grid, x2_grid, G, [1, 1], 'FaceColor', '#009900', 'EdgeColor', 'k'); % Green for target

% --- 5. Finalize the Plot ---
hold off; % Release the figure hold
grid on;
axis([x_1_min, x_1_max, x_2_min, x_2_max]); % Ensure axis limits are set correctly
% Set font sizes for axes and labels
set(gca, 'FontSize', axis_font_size);
xlabel('$\mathrm{x}$', 'FontSize', label_font_size);
ylabel('$\mathrm{v}$', 'FontSize', label_font_size);
xlim([-0.4,0.4])
ylim([-0.7,0.7])

% Add a matching color bar to show the value-to-color mapping
colormap(custom_cmap);
cb = colorbar;
ylabel(cb, 'g(s)', 'FontSize', label_font_size);
set(cb, 'FontSize', axis_font_size); % Set font size for colorbar ticks
caxis([0.05, 0.95]); % Set color bar limits to match the loop's range
% --- 6. Save the Figure -----
% Set the paper size to match the figure size to prevent labels from being cut off
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

% Save the figure to a PDF file with 500 DPI resolution
print(fig, 'g_IP', '-dpdf', '-r300');