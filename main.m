% AE 3030 Project 1 Thin Airfoil Theory
% Ashutosh Chakragiri and Hruday Shah
fprintf("Running AE 3030 Project 1:\n");

% Known variables
V_oo = 100;
p_oo = 101.3 * 1e3; % Pa
alpha = 15;
alpha_rad = deg2rad(15);
rho = 1.23; % kg / m^3

% Assumed variables
c = 1; % Chord length is 1 unit
n_vortices = 1000; % Number of vortices
n_streamlines = 100; % Number of streamlines
n_timesteps = 1600; % Number of timesteps
delta_t = 0.00003; % Size of each timestep

% Body code
vortex_spacing = c / n_vortices;

% Set up 100 evenly-spaced vortices between segments
vortices_x = linspace(0, c, n_vortices);
vortices_y = -vortices_x * sin(alpha_rad);

% Set up NACA 0000 airfoil with 15 deg. AoA
airfoil_x = 0:c/5:c*cos(alpha_rad);
airfoil_y = -airfoil_x * sin(alpha_rad);

fig = figure();
fig.WindowState = "maximized";

subplot(2, 2, [1, 3]);
plot(airfoil_x, airfoil_y, "black", "LineWidth", 2)
hold on;

% Plot vortex sheet strength
vortex_sheet_strength_x = 0:c/100:c;
vortex_sheet_strength_y = littleGamma(vortex_sheet_strength_x, alpha_rad, V_oo, c);

hold on;

% Plot strength of each vortex segment
vortex_segment_strength_y = zeros(1, length(vortices_x));

for i = 2:(length(vortices_x)-1)
    vortex_segment_strength_y(i) = bigGamma(vortices_x(i), alpha_rad, V_oo, c) - bigGamma(vortices_x(i - 1), alpha_rad, V_oo, c);
end

hold on;

% Calculate the lift coefficient
lift = rho * V_oo * sum(vortex_segment_strength_y);
c_l = lift / ((1 / 2) * rho * V_oo .^ 2 * c);
fprintf("- Coefficient of Lift: %d\n", c_l);

hold on;

% Initialize streamlines
streamline_pos_x = zeros(n_streamlines, n_timesteps);
streamline_pos_y = zeros(size(streamline_pos_x));
streamline_vel_mag = zeros(size(streamline_pos_x));
c_p = zeros(size(streamline_pos_x));

streamline_pos_x(:, 1) = -2 * c;
streamline_pos_y(:, 1) = linspace(-1.5, 1, n_streamlines);

% Calculate position of streamlines
for j = 2:n_timesteps
    u_vortex = 0;
    v_vortex = 0;
    for i = 1:n_vortices
        u_vortex = u_vortex + 1 / (2 * pi) * (vortex_segment_strength_y(i) * (streamline_pos_y(:, j - 1) - vortices_y(i))) ./ ((streamline_pos_x(:, j - 1) - vortices_x(i)) .^ 2 + (streamline_pos_y(:, j - 1) - vortices_y(i)) .^ 2);
        v_vortex = v_vortex - 1 / (2 * pi) * (vortex_segment_strength_y(i) * (streamline_pos_x(:, j - 1) - vortices_x(i))) ./ ((streamline_pos_x(:, j - 1) - vortices_x(i)) .^ 2 + (streamline_pos_y(:, j - 1) - vortices_y(i)) .^ 2);
    end
    streamline_vel_x = V_oo + u_vortex;
    streamline_vel_y = v_vortex;

    new_x = streamline_pos_x(:, j - 1) + streamline_vel_x * delta_t;
    new_y = streamline_pos_y(:, j - 1) + streamline_vel_y * delta_t;
    
    streamline_pos_x(:, j) = new_x;
    streamline_pos_y(:, j) = new_y;
    
    vel = sqrt(streamline_vel_x .^ 2 + streamline_vel_y .^ 2);
    streamline_vel_mag(:, j) = vel;
    c_p(:, j) = 1 - (vel / V_oo) .^ 2;
end

% Find stagnation point
[max_c_p, i] = max(c_p(:));
fprintf("- Coefficient of Pressure at Stag. Pt.: %d\n", max_c_p);
stag_pt = [streamline_pos_x(i), streamline_pos_y(i)];

plot(streamline_pos_x.', streamline_pos_y.', "blue");
hold on;
plot(stag_pt(1), stag_pt(2), "r*", "DisplayName", "Stagnation Point", "LineWidth", 1);
text(stag_pt(1) + 0.05, stag_pt(2), "Stagnation Point");

hold off;

% Set plot title and x-axis label
title("Airfoil with Streamlines vs. Distance");
ylabel("Vertical Distance (m)");
xlabel("Horizontal Distance (m)");
legend("Airfoil", "Streamline");
xlim([-2*c, 2*c]);

% Find c_p for streamline above and below stagnation point streamline
[row, col] = ind2sub(size(c_p), i);
c_p_row = 1;

% Plot streamline velocity and stagnation point
subplot(2, 2, 2);
plot(streamline_pos_x(row - c_p_row, :), streamline_vel_mag(row - c_p_row, :));
hold on;
plot(streamline_pos_x(row + c_p_row, :), streamline_vel_mag(row + c_p_row, :));

hold off;

% Set plot title and x-axis label
title("Velocity vs. Distance");
ylabel("Velocity (m/s)");
xlabel("Horizontal Distance (m)");
legend("V_{lower}", "V_{upper}");
xlim([-2*c, 2*c]);

% Get pressure distribution for one streamline above and below airfoil
p_dist_x = streamline_pos_x(row - c_p_row, :);
c_p_above = c_p(row - c_p_row, :);
c_p_below = c_p(row + c_p_row, :);

% Plot coefficient of pressure for one streamline above and below airfoil
subplot(2, 2, 4);
plot(streamline_pos_x(row - c_p_row, :), c_p_above);
hold on;
plot(streamline_pos_x(row + c_p_row, :), c_p_below);

h = gca;
set(h, "YDir", "reverse");

hold off;

% Set plot title and x-axis label
title("Coefficient of Pressure vs. Distance");
ylabel("c_p");
xlabel("Horizontal Distance (m)");
legend("Lower", "Upper");
xlim([-2*c, 2*c]);

fprintf(":---:\n");

% Functions
function integrated_circulation = bigGamma(s, alpha_rad, V_oo, c)
    integrated_circulation = 2 * alpha_rad * V_oo * (sqrt(s * (c - s)) - c * atan(sqrt((c - s) / s)));
end

function circ = littleGamma(s, alpha_rad, V_oo, c)
    circ = (2 .* alpha_rad .* V_oo .* sqrt((c - s) ./ s));
end
