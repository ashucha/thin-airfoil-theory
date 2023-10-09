%% Known variables
V_oo = 100;
alpha = 15;
alpha_rad = 15 * pi / 180;

%% Assumed variables
c = 100; % Chord length is 100 m

%% Body code

% Set up 9 evenly-spaced vortices between 10 segments
vortices_x = 5:10:c-5;
vortices_y = -(vortices_x - c / 2) * alpha_rad;

scatter(vortices_x, vortices_y, 50, "white");
hold on;
scatter(vortices_x, vortices_y, 25, "white");

hold on;

% Set up NACA 0000 airfoil with 15 deg. AoA
airfoil_x = [0, c];
airfoil_y = -(airfoil_x - c / 2) * alpha_rad;

plot(airfoil_x, airfoil_y, "white");
hold on;
plot(airfoil_x, zeros(2), "white", "Marker", ".");

% Plot vortex sheet strength
vortex_sheet_x = 0:2:c;
vortex_sheet_y = gamma(vortex_sheet_x);

plot(vortex_sheet_x, vortex_sheet_y, "cyan");

title("Project 1: Thin Airfoil Theory");
xlabel("Horizontal Distance (m)");

%% Functions
gamma = @(s) (2 .* alpha_rad .* V_oo .* sqrt((c - s) ./ s));

function strength = getVortexStrength(s1, s2)
    strength = integral(@gamma, s1, s2);
end
