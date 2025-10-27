%% MATLAB Code: Normalized Level Set Comparison (Final Tight Output with LaTeX Ticks)

clear
clc
close all
format long

% --- 0. Define File Paths and Common Parameters ---
ip_file_path = '/Users/saber/Desktop/Safety/IP/results/2025-10-09_13-55-18/AVR_gain_map.csv';
di_file_path = '/Users/saber/Desktop/Safety/DI/results/2025-10-08_23-15-25/AVR_gain_map.csv';
levels = 0.4:0.001:1;
num_levels = length(levels);

% --- 1. Processing Function (to avoid repeating code) ---
calculate_normalized_levels = @(g) computeLevelSet(g, levels, num_levels);

% --- 2. Process INVERTED PENDULUM (IP) Data (Error handling skipped for brevity) ---
try
    g_ip = readmatrix(ip_file_path);
    [y_ip, alpha_1_prop_ip] = calculate_normalized_levels(g_ip);
    success_ip = true;
catch
    warning('Could not load IP data from: %s', ip_file_path);
    y_ip = zeros(1, num_levels);
    success_ip = false;
end

% --- 3. Process DOUBLE INTEGRATOR (DI) Data (Error handling skipped for brevity) ---
try
    g_di = readmatrix(di_file_path);
    [y_di, alpha_1_prop_di] = calculate_normalized_levels(g_di);
    success_di = true;
catch
    warning('Could not load DI data from: %s', di_file_path);
    y_di = zeros(1, num_levels);
    success_di = false;
end

% --- 4. Prepare Plot Data ---
x_plot = fliplr(levels);
y_plot_ip = fliplr(y_ip);
y_plot_di = fliplr(y_di);

% --- 5. Generate and Customize Plot (Optimized for Tight Output) ---

% Set up figure and enable LaTeX rendering
fig = figure('Position', [100, 100, 800, 600], 'PaperPositionMode', 'auto'); 
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');
set(groot, 'DefaultAxesFontSize', 30);

ax = gca;
hold on;

% --- Plotting with LineWidth = 4 ---
if success_ip
    plot(x_plot, y_plot_ip, 'b-', 'LineWidth', 4, 'DisplayName', 'Inverted Pendulum');
end
if success_di
    plot(x_plot, y_plot_di, 'r--', 'LineWidth', 4, 'DisplayName', 'Double Integrator');
end
hold off;

% --- Title and Labels using LaTeX and increased FontSize ---
% title('Normalized Level Set Comparison', 'FontSize', 20);
xlabel('$\alpha$', 'FontSize', 30);
ylabel('$\frac{\mathcal{K}_{\alpha}}{\mathcal{K}_1}$', 'FontSize', 30, 'Interpreter', 'latex');

% *** CRUCIAL FIX: Ensure Tick Labels are Rendered in LaTeX ***
set(ax, 'TickLabelInterpreter', 'latex');

grid on;
legend('Location', 'northwest');

% --- Axis Configuration ---
xlim([0.4 1]); 
ylim([0.5,2])
set(gca, 'XDir', 'reverse'); 

% --- FIX: TIGHTENING THE PLOT ---
set(ax, 'units', 'centimeters');
ti = get(ax, 'TightInset');
op = get(ax, 'OuterPosition');
set(ax, 'units', 'normalized');
set(fig, 'PaperUnits', 'centimeters');
set(fig, 'PaperSize', [op(3)+ti(1)+ti(3), op(4)+ti(2)+ti(4)]); 
set(fig, 'PaperPosition', [0 0 op(3)+ti(1)+ti(3), op(4)+ti(2)+ti(4)]); 

% --- 6. Save Plot as PDF (Corrected for tight fit) ---
% Using '-painters' renderer for robust vector graphics and '-loose' for minimal margins.
print(fig, 'Normalized_Level_Set_Comparison', '-dpdf', '-r300', '-painters', '-loose'); 

disp('Plot generated and saved as "Normalized_Level_Set_Comparison.pdf" with all text (labels, title, and ticks) using LaTeX rendering.');

%% Local Function for Level Set Calculation
function [normalized_levels, alpha_1_proportion] = computeLevelSet(g, levels, num_levels)
    total_elements = numel(g);
    
    % Calculate the normalization factor
    alpha_1_size = nnz(g >= 1);
    alpha_1_proportion = alpha_1_size / total_elements;
    
    normalized_levels = zeros(1, num_levels);
    
    if alpha_1_proportion == 0
        % If normalization factor is zero, return the un-normalized proportion
        for i = 1:num_levels
            count = nnz(g >= levels(i));
            normalized_levels(i) = count / total_elements;
        end
    else
        % Normal case: calculate and normalize
        for i = 1:num_levels
            count = nnz(g >= levels(i));
            level_set_proportion = count / total_elements;
            normalized_levels(i) = level_set_proportion / alpha_1_proportion;
        end
    end
end