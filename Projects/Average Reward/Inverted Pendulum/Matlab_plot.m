clear all
clc
close all

% Set LaTeX interpreters
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultAxesTickLabelInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0, 'DefaultColorbarTickLabelInterpreter', 'latex');

% --- Updated Grid Settings ---
x_1_min = -0.5;
x_1_max =  0.5;
x_2_min = -1.0;
x_2_max =  1.0;
num_points_state_1 = 151;
num_points_state_2 = 151;

% --- Define Safe Box Region ---
box_x1_min = -0.3;
box_x1_max =  0.3;
box_x2_min = -0.6;
box_x2_max =  0.6;
safe_x = [box_x1_min, box_x1_max, box_x1_max, box_x1_min, box_x1_min];
safe_v = [box_x2_min, box_x2_min, box_x2_max, box_x2_max, box_x2_min];

dark_green = [0, 0.5, 0];

% --- Create Custom Red-to-Green Colormap ---
cmap_length = 256;
custom_cmap = zeros(cmap_length, 3);
for i = 1:cmap_length
    if i <= cmap_length/2
        custom_cmap(i, 1) = 1.0;
        custom_cmap(i, 2) = (i-1)/(cmap_length/2-1);
        custom_cmap(i, 3) = 0.0;
    else
        custom_cmap(i, 1) = 1.0 - (i-cmap_length/2)/(cmap_length/2);
        custom_cmap(i, 2) = 1.0;
        custom_cmap(i, 3) = 0.0;
    end
end

% --- Load Data ---
G_Average = readmatrix('results_primal_0.0/G.csv');

% --- Create Grid ---
x1_avg = linspace(x_1_min, x_1_max, num_points_state_1);
x2_avg = linspace(x_2_min, x_2_max, num_points_state_2);

% --- Apply Mask to Safe Box ---
rect_x_min = find(x1_avg >= box_x1_min, 1, 'first');
rect_x_max = find(x1_avg <= box_x1_max, 1, 'last');
rect_y_min = find(x2_avg >= box_x2_min, 1, 'first');
rect_y_max = find(x2_avg <= box_x2_max, 1, 'last');
mask = zeros(size(G_Average));
mask(rect_y_min:rect_y_max, rect_x_min:rect_x_max) = 1;
G_Average_masked = G_Average .* mask;
G_Average_masked(mask == 0) = NaN;

% --- Plotting ---
figure;
set(gcf, 'Position', [100, 100, 600, 500]);
hold on;

colormap(custom_cmap);
contour(x1_avg, x2_avg, G_Average_masked, 100);
caxis([0 0.9]);

% Colorbar
cb = colorbar;
ylabel(cb, 'Safety', 'FontSize', 14, 'Interpreter', 'latex');

% Highlight safe area where safety >= 0.99
high_safety_Avg = G_Average_masked >= 0.99;
[c_Avg, ~] = contour(x1_avg, x2_avg, high_safety_Avg, [0.5 0.5], '-k', 'LineWidth', 1.5);
if ~isempty(c_Avg)
    contour_data = c_Avg;
    level_idx = 1;
    while level_idx < size(contour_data, 2)
        n_points = contour_data(2, level_idx);
        x_coords = contour_data(1, level_idx+1:level_idx+n_points);
        y_coords = contour_data(2, level_idx+1:level_idx+n_points);
        fill(x_coords, y_coords, dark_green, 'EdgeColor', 'none');
        level_idx = level_idx + n_points + 1;
    end
end

% Draw the Safe Box
plot(safe_x, safe_v, 'k--', 'LineWidth', 2);

% Legend
rect_patch_Avg = patch('XData', [NaN], 'YData', [NaN], ...
    'FaceColor', dark_green, 'EdgeColor', 'none');
legend(rect_patch_Avg, {'Safe Set'}, 'Location', 'Northeast', ...
    'Interpreter', 'latex', 'FontSize', 12);

% Labels and Ticks
xlabel('$\mathrm{x}_1$', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('$\mathrm{x}_2$', 'FontSize', 20, 'FontWeight', 'bold');
title('Average Reward Safety Set', 'FontSize', 18, 'Interpreter', 'latex');

% Custom tick marks
x_tick_values = linspace(x_1_min, x_1_max, 5);
y_tick_values = linspace(x_2_min, x_2_max, 5);
set(gca, 'XTick', x_tick_values);
set(gca, 'YTick', y_tick_values);
set(gca, 'XTickLabel', arrayfun(@(x) ['$' num2str(x, '%.1f') '$'], x_tick_values, 'UniformOutput', false));
set(gca, 'YTickLabel', arrayfun(@(y) ['$' num2str(y, '%.1f') '$'], y_tick_values, 'UniformOutput', false));
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);

% Axis limits and layout
xlim([x_1_min, x_1_max]);
ylim([x_2_min, x_2_max]);
set(gcf, 'Color', 'white');

% % Export figure
% try
%     exportgraphics(gcf, 'Average_Reward_Safety_Boxed.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');
% catch
%     print('-dpdf', 'Average_Reward_Safety_Boxed.pdf');
% end
