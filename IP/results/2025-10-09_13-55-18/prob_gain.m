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
axis_font_size = 25;  % Font size for the tick labels (25)
label_font_size = 30; % Font size for the x/y labels (30)
% Define grid settings for the plot axes
x_1_min = -0.5; % Position (x) min
x_1_max = 0.5;  % Position (x) max
x_2_min = -1.0; % Velocity (v) min
x_2_max = 1.0;  % Velocity (v) max
% Define the rectangular safe set boundary for visualization
safe_x = [-0.3, 0.3, 0.3, -0.3, -0.3]; % Position coordinates
safe_v = [-0.6, -0.6, 0.6, 0.6, -0.6]; % Velocity coordinates
% --- 2. Data Loading and Preparation ---
% Load data from the specified CSV file
G = readmatrix('AVR_gain_map.csv');
% Create coordinate grid vectors corresponding to the dimensions of G
% x1_grid is Position (x) (Cols in G, now Y-axis)
% x2_grid is Velocity (v) (Rows in G, now X-axis)
x1_grid = linspace(x_1_min, x_1_max, size(G, 2)); % Position grid (x)
x2_grid = linspace(x_2_min, x_2_max, size(G, 1)); % Velocity grid (v)
% Create a mask to highlight the G=0 boundary area
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

% *** PLOT CHANGES: SWAP GRID VECTORS (x2_grid, x1_grid) AND TRANSPOSE DATA (G.') ***
% X-axis is now Velocity (x2_grid), Y-axis is now Position (x1_grid)
contourf(x2_grid, x1_grid, g_mask.', [0, 0], 'FaceColor', [1 0 0], 'EdgeColor', 'k'); % Grey for unsafe
contourf(x2_grid, x1_grid, g_mask.', [0.9999, 0.9999], 'FaceColor', 'g', 'EdgeColor', 'k'); % Green for target
% Plot the 999 intermediate contour levels with the color gradient
for i = 1:num_levels
    level = 0.0001 * i;
    line_color = custom_cmap(i, :); % Select the i-th color
    contour(x2_grid, x1_grid, G.', [level, level], 'LineWidth', 2.5, 'Color', line_color);
end
% *** PLOT CHANGE: SWAP SAFE BOUNDARY COORDINATES ***
% Plot takes (X, Y) -> (safe_v, safe_x)
plot(safe_v, safe_x, 'k--', 'LineWidth', 3);
contourf(x2_grid, x1_grid, G.', [1, 1], 'FaceColor', '#009900', 'EdgeColor', 'k'); % Green for target

% --- 5. Finalize the Plot ---
hold off; % Release the figure hold
grid on;

% *** AXIS LIMITS AND LABELS CHANGE: SWAP X and Y ***
% Ensure axis limits are set correctly
axis([x_2_min, x_2_max, x_1_min, x_1_max]); 
xlim([-0.7,0.7]) % X-limit for Velocity
ylim([-0.4,0.4]) % Y-limit for Position

set(gca, 'YTick', [-0.4, -0.2, 0.0, 0.2, 0.4]);

% Set font sizes and bold axes (from previous request)
set(gca, 'FontSize', axis_font_size, 'LineWidth', 3, 'Box', 'on'); 

 
% *** LABEL CHANGE: SWAP X and Y LABELS ***
xlabel('Velocity ($\mathrm{v}$)', 'FontSize', label_font_size); % X-axis is now Velocity
ylabel('Position ($\mathrm{x}$)', 'FontSize', label_font_size);  % Y-axis is now Position

% Add a matching color bar to show the value-to-color mapping
colormap(custom_cmap);
cb = colorbar;
ylabel(cb, 'g(s)', 'FontSize', label_font_size);
set(cb, 'FontSize', axis_font_size); % Set font size for colorbar ticks
caxis([0.05, 0.95]); % Set color bar limits to match the loop's range

% --- 6. Save the Figure (Fixed for Small File Size) -----
output_filename = 'g_IP.pdf'; % Output will be saved as g_IP.pdf

% Set the paper size for correct cropping (keeping your existing settings)
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

% Save the figure to a PDF file using the OpenGL renderer.
% The -opengl flag forces rasterization, which drastically reduces file size
% for plots with many vector objects (like your 9999 contours).
print(fig, output_filename, '-dpdf', '-r600', '-image'); 
fprintf('Figure saved as %s with rasterized content to reduce file size.\n', output_filename);