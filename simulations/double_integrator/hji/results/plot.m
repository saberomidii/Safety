%% Load CSV files
state = readmatrix('results/state_2d.csv');       % columns: x, v
V_matrix = readmatrix('value_function_gamma_201_girds.csv'); % size: [length(x1), length(x2)]

% L=3.605551275463989;
L = 0.0;
V_matrix = V_matrix + L;
% Extract coordinates
x = state(:,1);
v = state(:,2);

% Define grids according to how Julia saved the value function
x1 = unique(x);           % positions
x2 = unique(v);           % velocities

% Make sure V_matrix is in correct shape: x1 along rows, x2 along columns
if ~isequal(size(V_matrix), [length(x1), length(x2)])
    V_matrix = reshape(V_matrix, [length(x1), length(x2)]);
end

%% Plot: filled contour with box
box_x = 0; box_v = -3; box_width = 4; box_height = 6;

figure;
contourf(x1, x2, V_matrix, 10);  % transpose to match axes
colorbar;
xlabel('x'); ylabel('v');
title('Value Function Contour');
hold on;
rectangle('Position',[box_x, box_v, box_width, box_height], 'EdgeColor','k','LineWidth',2);
hold off;

%% Plot: zero level set with box
figure;
contour(x1, x2, V_matrix, [0. 0.], 'r', 'LineWidth', 2);  % zero level set
xlabel('x'); ylabel('v');
title('Zero Level Set of Value Function');
hold on;
rectangle('Position',[box_x, box_v, box_width, box_height], 'EdgeColor','k','LineWidth',2);
hold off;

%% Optional: 3D surface
figure;
surf(x1, x2, V_matrix);
shading interp;
colorbar;
xlabel('x'); ylabel('v'); zlabel('V');
title('Value Function Surface');


%% Extract zero-level set contour and compute area
C = contourc(x1, x2, V_matrix, [0.0 0.0]);  % get contour matrix for zero level set

total_area = 0;
k = 1;

while k < size(C,2)
    n = C(2,k);               % number of points in this contour segment
    x_contour = C(1, k+1:k+n);
    y_contour = C(2, k+1:k+n);
    
    % Area of polygon defined by this contour segment
    area_segment = polyarea(x_contour, y_contour);
    total_area = total_area + area_segment;
    
    k = k + n + 1;             % move to next contour segment
end

% Total area of state space
dx = x1(2)-x1(1);
dv = x2(2)-x2(1);
state_space_area = (max(x1)-min(x1))*(max(x2)-min(x2));

safe_set_percentage = total_area / state_space_area;





x1=linspace(-1,5,161);
x2=linspace(-5,5,161);