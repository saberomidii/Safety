clear all;
clc;
close all;
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
x_1_min = -1.0; % Position (x) min
x_1_max = 5.0;  % Position (x) max
x_2_min = -5.0; % Velocity (v) min
x_2_max = 5.0;  % Velocity (v) max
% Define the rectangular safe set boundary for visualization
safe_x = [0, 4, 4, 0, 0];   % Position coordinates
safe_v = [-3, -3, 3, 3, -3]; % Velocity coordinates
% --- 2. Data Loading and Preparation ---
% Load data from the specified CSV file
G = readmatrix('AVR_gain_map.csv');
% Create coordinate grid vectors corresponding to the dimensions of G
% x1_grid is typically the COLUMNS (X-axis in MATLAB's contour default) -> Now Y-axis (Position)
% x2_grid is typically the ROWS (Y-axis in MATLAB's contour default) -> Now X-axis (Velocity)
x1_grid = linspace(x_1_min, x_1_max, size(G, 2)); % Position grid (x)
x2_grid = linspace(x_2_min, x_2_max, size(G, 1)); % Velocity grid (v)
% Create a mask to highlight the G=0 boundary area
g_mask = G;
g_mask(1:32, 1:161) = -1;
g_mask(130:161, 1:161) = -1;
g_mask(32:130, 1:26) = -1;
g_mask(32:130, 136:161) = -1;
% --- 3. Create the Custom Colormap with Orange and Yellow ---
num_levels = 9999;
color_red    = [1, 0.65 0];
color_yellow = [1, 1, 0];
color_green  = [0.5, 1.0, 0.0];
num_steps_part1 = round(num_levels / 2);
num_steps_part2 = num_levels - num_steps_part1;
r1 = linspace(color_red(1), color_yellow(1), num_steps_part1);
g1 = linspace(color_red(2), color_yellow(2), num_steps_part1);
b1 = linspace(color_red(3), color_yellow(3), num_steps_part1);
part1 = [r1', g1', b1'];
r2 = linspace(color_yellow(1), color_green(1), num_steps_part2 + 1);
g2 = linspace(color_yellow(2), color_green(2), num_steps_part2 + 1);
b2 = linspace(color_yellow(3), color_green(3), num_steps_part2 + 1);
part2 = [r2(2:end)', g2(2:end)', b2(2:end)'];
custom_cmap = [part1; part2];
% --- 4. Generate the Plot ---
figure;
set(gcf, 'Position', [100, 100, 750, 550], 'Color', 'white');
hold on; 

% *** PLOT CHANGE: SWITCH X and Y GRID VECTORS, AND TRANSPOSE G ***
% The contour function is now called as: contourf(X_grid, Y_grid, Data_transposed, ...)
% Target: X-axis = Velocity (x2_grid), Y-axis = Position (x1_grid)
contourf(x2_grid, x1_grid, g_mask.', [0, 0], 'FaceColor', [1 0 0], 'EdgeColor', 'k');
contourf(x2_grid, x1_grid, g_mask.', [0.9999, 0.9999], 'FaceColor', 'g', 'EdgeColor', 'k');
for i = 1:num_levels
    level = 0.0001 * i;
    line_color = custom_cmap(i, :);
    % *** PLOT CHANGE: TRANSPOSE G (G.') ***
    contour(x2_grid, x1_grid, G.', [level, level], 'LineWidth', 2.5, 'Color', line_color);
end
% *** PLOT CHANGE: SWAP X and Y FOR SAFE BOUNDARY ***
% Plot takes (X, Y) -> (safe_v, safe_x)
plot(safe_v, safe_x, 'k--', 'LineWidth', 3);
contourf(x2_grid, x1_grid, G.', [1, 1], 'FaceColor', '#009900', 'EdgeColor', 'k');

% --- 5. Finalize the Plot ---
hold off; 
grid on;

% *** AXIS LIMITS AND LABELS CHANGE: SWAP X and Y ***
% OLD Limits (x1_min to x1_max for X, x2_min to x2_max for Y)
% NEW Limits (x2_min to x2_max for X, x1_min to x1_max for Y)
axis([x_2_min, x_2_max, x_1_min, x_1_max]);
xlim([-3.5, 3.5]); % Manually set limits for Velocity (v) on X-axis
ylim([-0.5, 4.5]); % Manually set limits for Position (x) on Y-axis

% Set font sizes and bold axes (from previous request)
set(gca, 'FontSize', axis_font_size, 'LineWidth', 3, 'Box', 'on'); 

% *** LABEL CHANGE: SWAP X and Y LABELS ***
xlabel('Velocity ($\mathrm{v}$)', 'FontSize', label_font_size); % X-axis is now Velocity
ylabel('Position ($\mathrm{x}$)', 'FontSize', label_font_size);  % Y-axis is now Position

% Add a matching color bar
colormap(custom_cmap);
cb = colorbar;
ylabel(cb, 'g(s)', 'FontSize', label_font_size); 
set(cb, 'FontSize', axis_font_size);
caxis([0.05, 0.95]);

% --- 6. Save the Figure (Fixed for Small File Size) -----
output_filename = 'g_di.pdf'; % Output will be saved as g_IP.pdf

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